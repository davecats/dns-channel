function mmap=storedstructure(fname,field,writable)

FSIZE=0;
for i=1:numel(field)/3
    if numel(field)/3 == 1; c=field{1}; else c=field{i,1}; end
    switch c
        case 'double'
            bs=8;
        case 'complex'
            bs=16;
        case 'uint8'
            bs=1;
    end
    if numel(field)/3 == 1; s=field{2}; else s=field{i,2}; end
    FSIZE=FSIZE+bs*prod(s);
end
CHUNKSIZE=100*1000
FSIZE

fh = fopen(fname,'r'); if fh>-1; fclose(fh); end
if fh==-1
    fh = fopen(fname,'w'); 
    for i=CHUNKSIZE:CHUNKSIZE:FSIZE;
        fwrite(fh,zeros(1,CHUNKSIZE),'uint8'); 
    end
    fwrite(fh,zeros(1,FSIZE-i),'uint8');
    fclose(fh);
end
mmap=memmapfile(fname, ...
                 'Format', field, ...
                 'Writable',writable);