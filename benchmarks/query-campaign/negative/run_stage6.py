"""Outlined-filter control: same original hot loop, no deferred Option boundary.
Includes an always-pass filter to quantify wrapper/code-generation overhead.
"""
import os,pathlib,subprocess
logs=pathlib.Path('campaign/logs');cpu=logs/'filter-negative-stage6-cpu.jsonl'
for label,index,n,kinds,traces in [
 ('gencode','campaign/index/gencode',10000,['pass','bloom8','fuse8'],[str(p) for p in sorted(pathlib.Path('campaign/negative-query').glob('*.trace'))]),
 ('gencode','campaign/index/gencode',100000,['pass','bloom8','fuse8'],[]),
 ('human','/scratch1/alevin-fry-ecosystem/tmp/sux-phastr/idx/human_r0152_t32',100000,['pass','bloom8','fuse8'],[]),
]:
 for kind in kinds:
  args=['campaign/negative-query/filter_compare-v4',index,str(n),'campaign/reads/rna100k_1.fq','0',f'campaign/negative-query/{label}-q0.keys',kind,*traces]
  name=f'filter-direct-{label}-n{n}-{kind}'
  os.link(cpu,logs/f'{name}-cpu.jsonl')
  with (logs/f'{name}.log').open('w') as out:subprocess.run(args,stdout=out,stderr=subprocess.STDOUT,check=True)
  print(label,n,kind,'complete',flush=True)
