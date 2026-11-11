"""Untimed, single-worker capture; output equality checked separately."""
import os,pathlib,subprocess
out=pathlib.Path('campaign/negative-query')
binary=out/'piscem-trace'
for dataset in ['sim100k','random100k','rna100k']:
 env=dict(os.environ,PISCEM_QUERY_TRACE=str(out/dataset))
 with (out/f'{dataset}-capture.log').open('w') as log:
  subprocess.run([str(binary),'map-bulk','--dict','sshash','--decoder','serial','-i','campaign/index/gencode','-1',str(out/f'{dataset}_1.fq'),'-2',str(out/f'{dataset}_2.fq'),'-o',str(out/f'{dataset}-captured'),'-t','1','-q'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 print(dataset,'captured',flush=True)
