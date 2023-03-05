% --- GraSPH Input Demo Matlab Script
% Script shows how the example hdf5 input dambreak.h5 can be generated using 
% Python and h5py. Each dataset is gzip compressed.h5. Note that Matlab stores
% arrays in column-major, same as Fortran.
% Note: you may get some warning about the type conversions clamping
% values (i.e., losing precision). This can be ignored.

filename = "dambreak.h5";

% Delete file before writing because Matlab hdf5 functions don't overwite
% files if they exist
delete(filename);

% constants
dim = 3; 
dx = 0.5;
rho = 1000;

% Preparing real particles' coordinate data
xp = 50; % number of particles in x direction
yp = 25; %                        y
zp = 50; %                        z
nreal = xp*yp*zp;
x = zeros(dim, nreal);
n = 0;
for i = 1:xp
    for j = 1:yp
        for k = 1:zp
            n = n + 1;
            x(1, n) = (i-0.5)*dx;
            x(2, n) = (j-0.5)*dx;
            x(3, n) = (k-0.5)*dx;
        end
    end
end

% write coordinate data
h5create(filename, "/real/x", [dim, nreal], "Datatype", "double", 'Deflate', 6, 'Chunksize', [dim, nreal]);
h5write(filename, "/real/x", x);

% write velocity data (zero)
h5create(filename, "/real/v", [dim, nreal], "Datatype", "double", 'Deflate', 6, 'Chunksize', [dim, nreal]);
h5write(filename, "/real/v", zeros(dim, nreal))

% write pressure data (zero)
h5create(filename, "/real/p", nreal, "Datatype", "double", 'Deflate', 6, 'Chunksize', nreal);
h5write(filename, "/real/p", zeros(1, nreal))

% write density data (rho)
h5create(filename, "/real/rho", nreal, "Datatype", "double", 'Deflate', 6, 'Chunksize', nreal);
h5write(filename, "/real/rho", linspace(rho, rho, nreal))

% write type data (one)
h5create(filename, "/real/type", nreal, "Datatype", "int32", 'Deflate', 6, 'Chunksize', nreal);
h5write(filename, "/real/type", ones(1, nreal))

% write global particle ID
h5create(filename, "/real/ind", nreal, "Datatype", "int32", 'Deflate', 6, 'Chunksize', nreal);
h5write(filename, "/real/ind", linspace(1, nreal, nreal))

% write number of real particles as attribute
h5writeatt(filename, '/real', 'n', int32(nreal));

% Preparing virtual particles' coordinate data
xp = 150; % number of particles in x direction
yp = 25; %                         y
zp = 100; %                         z
nlayer = 4;
nvirt = (xp+2*nlayer)*(yp+2*nlayer)*(zp+2*nlayer) - xp*yp*zp;
x = zeros(dim, nvirt);
itype = zeros(1, nvirt);
n = 0;
for i = 1-nlayer:xp+nlayer
    for j = 1-nlayer:yp+nlayer
        for k = 1-nlayer:zp+nlayer
            if i < 1 || j < 1 || k < 1 || i > xp || j > yp || k > zp
                n = n + 1;
                x(1, n) = (i-0.5)*dx;
                x(2, n) = (j-0.5)*dx;
                x(3, n) = (k-0.5)*dx;
                if (i < 1 && j < 1) || (i < 1 && j > yp) || (i > xp && j < 1) || (i > xp && j > yp)
                    itype(n) = -5;
                elseif j < 1 || j > yp
                    itype(n) = -4;
                elseif i < 1 || i > xp
                    itype(n) = -3;
                elseif k > zp
                    itype(n) = -2;
                elseif k < 1
                    itype(n) = -1;
                end
            end
        end
    end
end

% write coordinate data
h5create(filename, "/virt/x", [dim, nvirt], "Datatype", "double", 'Deflate', 6, 'Chunksize', [dim, nvirt]);
h5write(filename, "/virt/x", x);

% write velocity data (zero)
h5create(filename, "/virt/v", [dim, nvirt], "Datatype", "double", 'Deflate', 6, 'Chunksize', [dim, nvirt]);
h5write(filename, "/virt/v", zeros(dim, nvirt))

% write pressure data (zero)
h5create(filename, "/virt/p", nvirt, "Datatype", "double", 'Deflate', 6, 'Chunksize', nvirt);
h5write(filename, "/virt/p", zeros(1, nvirt))

% write density data (rho)
h5create(filename, "/virt/rho", nvirt, "Datatype", "double", 'Deflate', 6, 'Chunksize', nvirt);
h5write(filename, "/virt/rho", linspace(rho, rho, nvirt))

% write type data (one)
h5create(filename, "/virt/type", nvirt, "Datatype", "int32", 'Deflate', 6, 'Chunksize', nvirt);
h5write(filename, "/virt/type", itype)

% write global particle ID
h5create(filename, "/virt/ind", nvirt, "Datatype", "int32", 'Deflate', 6, 'Chunksize', nvirt);
h5write(filename, "/virt/ind", linspace(1, nvirt, nvirt))

% write number of real particles as attribute
h5writeatt(filename, '/virt', 'n', int32(nvirt));