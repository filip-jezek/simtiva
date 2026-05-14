# Generate PDF from HTML

Convert an HTML file to PDF using weasyprint.

## Usage

```
/pdf [filename.html]
```

If no filename is provided, the most recently modified `.html` file in the current directory is used.

## Steps

1. Determine the input HTML file:
   - If `$ARGUMENTS` is provided and non-empty, use it as the path
   - Otherwise, run: `find . -maxdepth 1 -name "*.html" -printf "%T@ %p\n" | sort -rn | head -1 | awk '{print $2}'`

2. Derive the output PDF path by replacing the `.html` extension with `.pdf`

3. Run weasyprint:
   ```bash
   python3 -m weasyprint <input.html> <output.pdf>
   ```

4. Report the output PDF path and file size:
   ```bash
   ls -lh <output.pdf>
   ```

If weasyprint is not installed, install it first:
```bash
pip install weasyprint
```
