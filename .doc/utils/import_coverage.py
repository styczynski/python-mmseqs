def import_coverage_report():
    import glob
    import os
    from shutil import copyfile

    from pathlib import Path
    Path("./_static").mkdir(parents=True, exist_ok=True)
    Path("./_sphinx_resources").mkdir(parents=True, exist_ok=True)

    replacements = dict()
    html_contents_dict = dict()
    html_files = set()

    for html_file in glob.glob('../test_reports/coverage_latest_html/*.html', recursive=True):
        html_files.add(html_file)
        with open(html_file, 'r') as in_html:
            html_contents = in_html.read()
            html_contents = html_contents.replace('<body>', '').replace('</body>', '').replace('<head>', '').replace('</head>', '')
            html_contents = html_contents.replace('<!DOCTYPE html>', '').replace('<html>', '').replace('</html>', '')

        html_contents = f"""
<style>
    #keyboard_icon {{
        width: 55px;
        margin-left: 170px !important;
        position: relative;
        top: 27px;
    }}
</style>
{html_contents}"""

        name = 'coverage_'+os.path.basename(html_file).replace('.html', '')
        replacements[os.path.basename(html_file)] = f'{name}.html'
        out_rst_path = f'_sphinx_resources/{name}.rst'
        out_html_path = f'_sphinx_resources/{name}.html'

        rst_contents = f"""Coverage {name}
============================

.. raw:: html
   :file: {name}.html
"""

        with open(out_rst_path, 'w') as out_rst:
            out_rst.write(rst_contents)
        html_contents_dict[out_html_path] = html_contents

    for asset_file in glob.glob('../test_reports/coverage_latest_html/*', recursive=True):
        if asset_file not in html_files:
            if os.path.isfile(asset_file):
                copyfile(asset_file, f'./_static/coverage_{os.path.basename(asset_file)}')
                replacements[os.path.basename(asset_file)] = f'../_static/coverage_{os.path.basename(asset_file)}'

    for out_html_path in html_contents_dict.keys():
        content = html_contents_dict[out_html_path]
        with open(out_html_path, 'w') as out_html:
            for to_repalce in replacements.keys():
                content = content.replace(to_repalce, replacements[to_repalce])
            out_html.write(content)
