use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use syn::visit::{self, Visit};
use syn::{ItemUse, UseTree};

const LAYERS: [&str; 5] = ["domain", "adapters", "engines", "application", "cli_compat"];

fn rust_files(root: &Path) -> Vec<PathBuf> {
    let mut pending = vec![root.to_path_buf()];
    let mut files = Vec::new();
    while let Some(path) = pending.pop() {
        for entry in fs::read_dir(&path).expect("architecture directory is readable") {
            let path = entry.expect("valid directory entry").path();
            if path.is_dir() {
                pending.push(path);
            } else if path.extension().is_some_and(|extension| extension == "rs") {
                files.push(path);
            }
        }
    }
    files.sort();
    files
}

fn module_path(source_root: &Path, path: &Path) -> Vec<String> {
    let relative = path
        .strip_prefix(source_root)
        .expect("source is below crate root");
    let mut components = relative
        .components()
        .map(|component| component.as_os_str().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let file = components.pop().expect("Rust source has a file name");
    let stem = file.strip_suffix(".rs").expect("Rust source suffix");
    if stem != "mod" {
        components.push(stem.to_string());
    }
    components
}

fn flatten_use(tree: &UseTree, prefix: &mut Vec<String>, output: &mut Vec<Vec<String>>) {
    match tree {
        UseTree::Path(path) => {
            prefix.push(path.ident.to_string());
            flatten_use(&path.tree, prefix, output);
            prefix.pop();
        }
        UseTree::Name(name) => {
            prefix.push(name.ident.to_string());
            output.push(prefix.clone());
            prefix.pop();
        }
        UseTree::Rename(rename) => {
            prefix.push(rename.ident.to_string());
            output.push(prefix.clone());
            prefix.pop();
        }
        UseTree::Glob(_) => output.push(prefix.clone()),
        UseTree::Group(group) => {
            for item in &group.items {
                flatten_use(item, prefix, output);
            }
        }
    }
}

struct DependencyVisitor {
    current: Vec<String>,
    references: Vec<Vec<String>>,
}

impl<'ast> Visit<'ast> for DependencyVisitor {
    fn visit_item_use(&mut self, node: &'ast ItemUse) {
        flatten_use(&node.tree, &mut Vec::new(), &mut self.references);
        visit::visit_item_use(self, node);
    }

    fn visit_path(&mut self, node: &'ast syn::Path) {
        self.references.push(
            node.segments
                .iter()
                .map(|segment| segment.ident.to_string())
                .collect(),
        );
        visit::visit_path(self, node);
    }
}

fn canonical_reference(current: &[String], raw: &[String]) -> Option<String> {
    if raw.is_empty() {
        return None;
    }
    let mut resolved = current.to_vec();
    let mut index = 0;
    match raw[0].as_str() {
        "crate" => {
            resolved.clear();
            index = 1;
        }
        "self" => index = 1,
        "super" => {
            while raw.get(index).is_some_and(|segment| segment == "super") {
                resolved.pop();
                index += 1;
            }
        }
        first if LAYERS.contains(&first) => resolved.clear(),
        _ => return None,
    }
    resolved.extend(raw[index..].iter().cloned());
    if resolved.len() < 2 || !LAYERS.contains(&resolved[0].as_str()) {
        return None;
    }
    Some(format!("{}::{}", resolved[0], resolved[1]))
}

fn owner(module: &[String]) -> String {
    format!(
        "{}::{}",
        module[0],
        module.get(1).map(String::as_str).unwrap_or("__root")
    )
}

fn allowed_dependency(from: &str, to: &str) -> bool {
    let from_layer = from.split("::").next().unwrap();
    let to_layer = to.split("::").next().unwrap();
    match from_layer {
        "domain" => to_layer == "domain",
        "adapters" => matches!(to_layer, "domain" | "adapters"),
        "engines" => matches!(to_layer, "domain" | "adapters" | "engines"),
        "application" => true,
        "cli_compat" => matches!(to_layer, "domain" | "cli_compat"),
        _ => false,
    }
}

fn dependency_graph() -> BTreeMap<String, BTreeSet<String>> {
    let source_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("rust/src");
    let mut graph = BTreeMap::<String, BTreeSet<String>>::new();
    for layer in LAYERS {
        for path in rust_files(&source_root.join(layer)) {
            let module = module_path(&source_root, &path);
            let source_owner = owner(&module);
            let source = fs::read_to_string(&path).expect("Rust source is readable");
            let syntax = syn::parse_file(&source)
                .unwrap_or_else(|error| panic!("failed to parse {}: {error}", path.display()));
            let mut visitor = DependencyVisitor {
                current: module,
                references: Vec::new(),
            };
            visitor.visit_file(&syntax);
            let edges = graph.entry(source_owner.clone()).or_default();
            for raw in visitor.references {
                if let Some(target) = canonical_reference(&visitor.current, &raw)
                    && target != source_owner
                {
                    assert!(
                        allowed_dependency(&source_owner, &target),
                        "forbidden architectural dependency {source_owner} -> {target} in {}",
                        path.display()
                    );
                    edges.insert(target);
                }
            }
        }
    }
    graph
}

fn visit_node(
    node: &str,
    graph: &BTreeMap<String, BTreeSet<String>>,
    visiting: &mut Vec<String>,
    complete: &mut BTreeSet<String>,
) {
    if complete.contains(node) {
        return;
    }
    if let Some(index) = visiting.iter().position(|candidate| candidate == node) {
        let mut cycle = visiting[index..].to_vec();
        cycle.push(node.to_string());
        panic!("architectural dependency cycle: {}", cycle.join(" -> "));
    }
    visiting.push(node.to_string());
    if let Some(edges) = graph.get(node) {
        for target in edges {
            visit_node(target, graph, visiting, complete);
        }
    }
    visiting.pop();
    complete.insert(node.to_string());
}

#[test]
fn canonical_layer_graph_is_acyclic_and_points_inward() {
    let graph = dependency_graph();
    let mut complete = BTreeSet::new();
    for node in graph.keys() {
        visit_node(node, &graph, &mut Vec::new(), &mut complete);
    }
}

#[test]
fn crate_root_exposes_only_run() {
    let source = fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join("rust/src/lib.rs"))
        .expect("crate root is readable");
    let syntax = syn::parse_file(&source).expect("crate root parses");
    let public = syntax
        .items
        .iter()
        .filter_map(|item| match item {
            syn::Item::Fn(function) if matches!(function.vis, syn::Visibility::Public(_)) => {
                Some(function.sig.ident.to_string())
            }
            syn::Item::Mod(module) if matches!(module.vis, syn::Visibility::Public(_)) => {
                Some(module.ident.to_string())
            }
            syn::Item::Use(item) if matches!(item.vis, syn::Visibility::Public(_)) => {
                Some("public use".to_string())
            }
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(public, ["run"]);
}
