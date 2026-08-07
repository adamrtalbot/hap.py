//! Pure genotype projections used by preprocessing orchestration.

pub(super) fn project_split_genotype(gt: &str, target: usize) -> String {
    let separator = if gt.contains('|') { '|' } else { '/' };
    let alleles = gt.split(['/', '|']).collect::<Vec<_>>();
    if alleles.contains(&".") {
        return gt.to_string();
    }
    let called = alleles
        .iter()
        .filter(|allele| allele.parse::<usize>().ok() == Some(target))
        .count();
    match called {
        0 => vec!["0"; alleles.len()].join(&separator.to_string()),
        count if count == alleles.len() => vec!["1"; alleles.len()].join(&separator.to_string()),
        _ => ["0", "1"].join(&separator.to_string()),
    }
}

pub(super) fn project_split_ad(ad: &str, target: usize) -> String {
    if ad == "." || ad.is_empty() {
        return ad.to_string();
    }
    let values = ad.split(',').collect::<Vec<_>>();
    if values.len() <= 2 {
        return ad.to_string();
    }
    format!(
        "{},{}",
        values.first().copied().unwrap_or("0"),
        values.get(target).copied().unwrap_or("0")
    )
}

pub(super) fn remap_gt(gt: &str, mapping: &[usize]) -> String {
    let separator = if gt.contains('|') { "|" } else { "/" };
    gt.split(['/', '|'])
        .map(|token| {
            token.parse::<usize>().ok().map_or_else(
                || token.to_string(),
                |allele| mapping.get(allele).copied().unwrap_or(0).to_string(),
            )
        })
        .collect::<Vec<_>>()
        .join(separator)
}

pub(super) fn bcf_encoded_gt(gt: &str) -> String {
    let mut phased_next = false;
    let mut values = Vec::new();
    let mut token = String::new();
    let flush = |token: &mut String, phased: bool, values: &mut Vec<String>| {
        if token.is_empty() {
            return;
        }
        let encoded = if token == "." {
            0
        } else {
            token
                .parse::<i32>()
                .map(|allele| ((allele + 1) << 1) | i32::from(phased))
                .unwrap_or(0)
        };
        values.push(encoded.to_string());
        token.clear();
    };
    for character in gt.chars() {
        match character {
            '/' | '|' => {
                flush(&mut token, phased_next, &mut values);
                phased_next = character == '|';
            }
            _ => token.push(character),
        }
    }
    flush(&mut token, phased_next, &mut values);
    values.join(",")
}

/// Expand a haploid GT token to its diploid legacy form.
pub(super) fn expand_haploid_gt(gt: &str, sex_chromosome: bool) -> String {
    if let Some(separator) = gt.chars().find(|separator| matches!(separator, '/' | '|')) {
        let alleles = gt.split(separator).collect::<Vec<_>>();
        if alleles.len() == 2
            && alleles[0] == "."
            && alleles[1].parse::<u32>().is_ok_and(|allele| allele > 0)
        {
            return format!("0{separator}{}", alleles[1]);
        }
        return gt.to_string();
    }
    match gt {
        "." => "./.".to_string(),
        "0" => "0/0".to_string(),
        other => match other.parse::<u32>() {
            Ok(n) if n > 0 && sex_chromosome => format!("{n}/{n}"),
            Ok(n) if n > 0 => format!("0/{n}"),
            _ => gt.to_string(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::{
        bcf_encoded_gt, expand_haploid_gt, project_split_ad, project_split_genotype, remap_gt,
    };

    #[test]
    fn projects_genotypes_without_vcf_or_filesystem_state() {
        assert_eq!(project_split_genotype("2/1", 2), "0/1");
        assert_eq!(project_split_ad("10,3,7", 2), "10,7");
        assert_eq!(remap_gt("2|1", &[0, 2, 1]), "1|2");
        assert_eq!(expand_haploid_gt("2", false), "0/2");
        assert_eq!(bcf_encoded_gt("0|1"), "2,5");
    }
}
