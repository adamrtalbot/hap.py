// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";
import catppuccin from "@catppuccin/starlight";

const site = process.env.SITE_URL || "https://adamrtalbot.github.io";
const base = (process.env.BASE_PATH || "/hap.py").replace(/\/$/, "");

export default defineConfig({
  site,
  base,
  integrations: [
    starlight({
      title: "hap-rs",
      description: "Haplotype-aware variant benchmarking in one native executable.",
      logo: {
        src: "./src/assets/hap-rs-mark.svg",
        alt: "hap-rs",
      },
      favicon: "/favicon.svg",
      head: [
        {
          tag: "meta",
          attrs: {
            property: "og:image",
            content: `${site}${base}/og.png`,
          },
        },
      ],
      social: [
        {
          icon: "github",
          label: "GitHub",
          href: "https://github.com/adamrtalbot/hap.py",
        },
      ],
      editLink: {
        baseUrl: "https://github.com/adamrtalbot/hap.py/edit/master/docs/src/content/docs/",
      },
      sidebar: [
        {
          label: "Getting Started",
          items: [
            { label: "Introduction", slug: "getting-started/introduction" },
            { label: "Installation", slug: "getting-started/installation" },
            { label: "Quick Start", slug: "getting-started/quick-start" },
            { label: "Core Concepts", slug: "getting-started/concepts" },
          ],
        },
        {
          label: "Tools",
          items: [
            { label: "Germline", slug: "tools/germline" },
            { label: "Somatic", slug: "tools/somatic" },
            { label: "Preprocess", slug: "tools/pre" },
            { label: "Feature Extraction", slug: "tools/ftx" },
            { label: "Quantify", slug: "tools/quantify" },
            { label: "Validate", slug: "tools/validate" },
          ],
        },
        {
          label: "Reference",
          items: [
            { label: "Commands & Engines", slug: "reference/commands-engines" },
            { label: "Inputs & Outputs", slug: "reference/inputs-outputs" },
            { label: "Metrics", slug: "reference/metrics" },
            { label: "Troubleshooting", slug: "reference/troubleshooting" },
          ],
        },
        {
          label: "Project",
          items: [
            { label: "Verification", slug: "project/verification" },
            { label: "Contributing", slug: "project/contributing" },
            { label: "License & Attribution", slug: "project/license" },
          ],
        },
      ],
      customCss: ["./src/styles/custom.css"],
      components: {
        Footer: "./src/components/Footer.astro",
      },
      plugins: [
        catppuccin({
          dark: { flavor: "mocha", accent: "teal" },
          light: { flavor: "latte", accent: "teal" },
        }),
      ],
    }),
  ],
});
