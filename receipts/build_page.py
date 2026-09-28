#!/usr/bin/env python3
"""build_page.py -- inline rh_calculator_embed.json and rh_math.js into index.template.html -> ../index.html (self-contained).
Run rh_calculator_twin.py first (writes the embed), then `node check_page_math.js` (headless self-check). Prints the sha256 of the built page."""
import hashlib, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
tpl = open(os.path.join(HERE, 'index.template.html'), encoding='utf-8').read()
embed = open(os.path.join(HERE, 'rh_calculator_embed.json'), encoding='utf-8').read()
mathjs = open(os.path.join(HERE, 'rh_math.js'), encoding='utf-8').read()
json.loads(embed)
assert tpl.count('__EMBED_JSON__') == 1 and tpl.count('__RH_MATH_JS__') == 1
assert '</script>' not in mathjs and '</script>' not in embed
html = tpl.replace('__EMBED_JSON__', embed).replace('__RH_MATH_JS__', mathjs)
out = os.path.join(HERE, '..', 'index.html')
open(out, 'w', encoding='utf-8').write(html)
print(out, len(html), 'bytes  sha256', hashlib.sha256(html.encode()).hexdigest())
