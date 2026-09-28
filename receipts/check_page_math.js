#!/usr/bin/env node
// check_page_math.js — run the page's own self-check headlessly: the same rh_math.js the browser executes,
// against the same embedded reference values the Python twin wrote. Usage: node check_page_math.js
'use strict';
const fs = require('fs'), path = require('path');
require(path.join(__dirname, 'rh_math.js'));
const EMB = JSON.parse(fs.readFileSync(path.join(__dirname, 'rh_calculator_embed.json'), 'utf8'));
RH.setZeros(EMB.zeros);
const t0 = Date.now();
const r = RH.selfCheck(EMB.ref);
for (const it of r.items) if (!it.pass) console.log('✗', it.name, '\n   page', JSON.stringify(it.got), '\n   twin', JSON.stringify(it.want));
console.log(`self-check ${r.pass}/${r.total} against the Python twin (${((Date.now()-t0)/1000).toFixed(1)} s)`);
process.exit(r.pass === r.total ? 0 : 1);
