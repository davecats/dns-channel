function mmap=storedstructure(fname,field,writable,create)

FSIZE=0;
for i=1:length(field{:,1})
    switch field{i,end}
        case 'double'
            bs=8;
        case 'complex'
            bs=16;
        case 'uint8'
            bs=1;
    end
    FSIZE=FSIZE+bs*prod(field{i,2});
end

fh = fopen(fname,'w'); 
for iy=0:dns.ny
    fwrite(fh,zeros(2,2*dns.nz+1,2*dns.nx+1),'double'); 
end
fclose(fh);
PSDimage=memmapfile(PSD.fname, ...
                 'Format', {'double' [2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'PSD'; ...
                            'double' [2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'CORR'}, ...
                 'Writable',true);