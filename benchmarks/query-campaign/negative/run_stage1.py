"""Run stage 1 sequentially under the quiet runner; retain common telemetry."""
import os,pathlib,subprocess
logs=pathlib.Path('campaign/logs')
cpu=logs/'filter-negative-stage1-cpu.jsonl'
for label,index,traces in [
 ('gencode','campaign/index/gencode',[str(p) for p in sorted(pathlib.Path('campaign/negative-query').glob('*.trace'))]),
 ('human','/scratch1/alevin-fry-ecosystem/tmp/sux-phastr/idx/human_r0152_t32',[]),
]:
 name=f'filter-negative-pipeline-{label}'
 os.link(cpu,logs/f'{name}-cpu.jsonl')
 with (logs/f'{name}.log').open('w') as out:
  subprocess.run(['campaign/negative-query/negative_pipeline',index,'10000','campaign/reads/rna100k_1.fq',*traces],stdout=out,stderr=subprocess.STDOUT,check=True)
 print(label,'complete',flush=True)
