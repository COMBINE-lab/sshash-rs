"""Full Gencode filters; one construction per process, rotated family order."""
import os,pathlib,subprocess
logs=pathlib.Path('campaign/logs');cpu=logs/'filter-negative-stage3-cpu.jsonl'
for q,kinds in [(16,['fuse8','zor4','bloom8','bitmap','zor5','bloom16']),(18,['bloom16','zor5','bitmap','bloom8','zor4','fuse8'])]:
 keys=f'campaign/negative-query/gencode-q{q}.keys'
 args=['campaign/negative-query/filter_compare','campaign/index/gencode','10000','campaign/reads/rna100k_1.fq',str(q),keys]
 if not pathlib.Path(keys).exists():
  with (logs/f'extract-gencode-q{q}.log').open('w') as out:subprocess.run([*args,'extract'],stdout=out,stderr=subprocess.STDOUT,check=True)
 for kind in kinds:
  name=f'filter-real-gencode-q{q}-{kind}'
  os.link(cpu,logs/f'{name}-cpu.jsonl')
  with (logs/f'{name}.log').open('w') as out:subprocess.run([*args,kind,*map(str,sorted(pathlib.Path('campaign/negative-query').glob('*.trace')))],stdout=out,stderr=subprocess.STDOUT,check=True)
  print(q,kind,'complete',flush=True)
