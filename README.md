# Allyson DeMerlis' Online Lab Notebook

Source for [ademerlis.github.io](https://ademerlis.github.io), built with [Quarto](https://quarto.org). My research website is [allysondemerlis.com](https://allysondemerlis.com).

## Adding a post

1. Run `./new_post.sh` (or copy any file in `posts/`) to create `posts/YYYY-MM-DD-short-title.md`.
2. Fill in the title and categories at the top, then write the post in Markdown.
3. Commit and push to `main`. GitHub Actions rebuilds and publishes the site in about two minutes.

To preview locally before pushing, run `quarto preview` in this folder.

## Layout

- `posts/`: one Markdown file per notebook entry
- `index.qmd`, `notebook.qmd`, `projects.qmd`, `about.qmd`: main pages
- `_quarto.yml`: site settings and navigation
- `styles.scss`: colors and fonts
- `.github/workflows/publish.yml`: builds and deploys the site on every push to `main`
