"""Net pre-MPHF streaming tests plus larger-query-set cache sensitivity.
Same v3 binary and scalar controls throughout. No production API is replaced.
"""
import os,pathlib,subprocess
logs=pathlib.Path('campaign/logs');cpu=logs/'filter-negative-stage5-cpu.jsonl'
for label,index,n,cases,traces in [
 ('gencode','campaign/index/gencode',10000,[(0,'bloom8'),(0,'fuse8')],[str(p) for p in sorted(pathlib.Path('campaign/negative-query').glob('*.trace'))]),
 ('gencode','campaign/index/gencode',100000,[(0,'fuse8'),(0,'bloom8'),(16,'fuse8'),(16,'bitmap'),(18,'fuse8'),(18,'bloom8')],[]),
 ('human','/scratch1/alevin-fry-ecosystem/tmp/sux-phastr/idx/human_r0152_t32',100000,[(0,'bloom8'),(0,'fuse8'),(16,'bitmap'),(18,'fuse8')],[]),
]:
 for q,kind in cases:
  keys=f'campaign/negative-query/{label}-q{q}.keys'
  args=['campaign/negative-query/filter_compare-v3',index,str(n),'campaign/reads/rna100k_1.fq',str(q),keys,kind,*traces]
  name=f'filter-net-{label}-n{n}-q{q}-{kind}'
  os.link(cpu,logs/f'{name}-cpu.jsonl')
  with (logs/f'{name}.log').open('w') as out:subprocess.run(args,stdout=out,stderr=subprocess.STDOUT,check=True)
  print(label,n,q,kind,'complete',flush=True)
