"""Selected-minimizer pre-MPHF filter comparison and human q-mer finalists.
This follows, rather than overlaps, the initial full Gencode filter campaign.
"""
import os,pathlib,subprocess
logs=pathlib.Path('campaign/logs');cpu=logs/'filter-negative-stage4-cpu.jsonl'
for label,index,qs,traces in [
 ('gencode','campaign/index/gencode',[0,16,18],[str(p) for p in sorted(pathlib.Path('campaign/negative-query').glob('*.trace'))]),
 ('human','/scratch1/alevin-fry-ecosystem/tmp/sux-phastr/idx/human_r0152_t32',[0,16,18],[]),
]:
 for q in qs:
  keys=f'campaign/negative-query/{label}-q{q}.keys'
  args=['campaign/negative-query/filter_compare-v2',index,'10000','campaign/reads/rna100k_1.fq',str(q),keys]
  if not pathlib.Path(keys).exists():
   with (logs/f'extract-{label}-q{q}.log').open('w') as out:subprocess.run([*args,'extract'],stdout=out,stderr=subprocess.STDOUT,check=True)
  kinds=['mphf7','fuse8','bloom8','zor4','bloom16','zor5'] if q==0 else ['fuse8','bloom8','bitmap']
  for kind in kinds:
   name=f'filter-confirm-{label}-q{q}-{kind}' if label=='gencode' and q else f'filter-real-{label}-q{q}-{kind}'
   os.link(cpu,logs/f'{name}-cpu.jsonl')
   with (logs/f'{name}.log').open('w') as out:subprocess.run([*args,kind,*traces],stdout=out,stderr=subprocess.STDOUT,check=True)
   print(label,q,kind,'complete',flush=True)
