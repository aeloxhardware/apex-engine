import json, sys
rows = [json.loads(l) for l in open(r'D:\apex_multibranch\results.jsonl')]
ref = {r['model']: r for r in rows if r['variant'] == 'ref' and 'status' not in r}
models = ['real_overhang.stl','dino.stl','benchy.stl','death.stl','fig.stl','l_shape.stl','l_triple.stl','deer.stl']
variants = []
for r in rows:
    if r['variant'] not in variants and r['variant'] != 'ref':
        variants.append(r['variant'])
latest = {}
for r in rows:
    latest[(r['model'], r['variant'])] = r
def cell(m, v):
    r = latest.get((m, v))
    if r is None: return '   -   '
    if 'status' in r: return r['status'][:9]
    b = ref[m]
    return '%5.2fx|%+.3f' % (r['unsupported_fraction'] / max(b['unsupported_fraction'], 1e-9),
                             (r.get('collision_fraction') or 0) - (b.get('collision_fraction') or 0))
print('format : ratio matiere sans appui (>1 = pire) | delta fraction de collisions')
print('%-18s' % 'variante' + ''.join('%-16s' % m[:14] for m in models))
print('%-18s' % 'ref(abs %/coll)' + ''.join('%-16s' % ('%.2f%%/%.3f' % (100 * ref[m]['unsupported_fraction'], ref[m].get('collision_fraction') or 0)) for m in models))
for v in variants:
    print('%-18s' % v + ''.join('%-16s' % cell(m, v) for m in models))


