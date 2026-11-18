#!/usr/bin/env python3
"""Record installed pyfastx semantics after a missing query; no baseline patch."""
import json
from pathlib import Path
import sys
import pyfastx
root=Path(sys.argv[1]);results=[]
for fastq in (False,True):
    d=root/'smoke-missing-failed'/str(int(fastq))
    for form in ('plain','gzip'):
        kind=pyfastx.Fastq if fastq else pyfastx.Fasta
        index=kind(str(d/('source.gz' if form=='gzip' else 'source')),index_file=str(d/f'{form}.fxi'),**({} if fastq else {'uppercase':False}))
        values=[]
        for key in ('last','missing','last','id:1/2'):
            try:
                r=index[key];values.append({'id':key,'sequence':r.seq,'quality':r.qual if fastq else ''})
            except KeyError as e:values.append({'id':key,'error':str(e)})
        results.append({'fastq':fastq,'form':form,'queries':values})
print(json.dumps({'version':pyfastx.__version__,'results':results},indent=2))
