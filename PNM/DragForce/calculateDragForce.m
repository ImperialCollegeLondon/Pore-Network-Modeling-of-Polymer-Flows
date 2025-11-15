function [PN_force, id_complete] = calculateDragForce(CasePath, DEMInfo, P_pore, throatConn, dx)

networkPath = [CasePath, 'pore_network/'];
spInterface = readSolidPoreInterface([CasePath, 'dual_network_interface_packing.csv'], dx);
poreInfo = readCenterAndSize([CasePath, 'pore_center_packing.csv'], dx);
particleInfo = readCenterAndSize([CasePath, 'solid_center_packing.csv'], dx);
index_list = dlmread([CasePath, 'index_list.csv'], ',', 1,0); index_list = index_list(2:end, 2:3);
completeParticleId = dlmread([CasePath, 'completeSpheres.dat'], ',');
for i = 1:length(completeParticleId)
    id_complete(i) = find(index_list(:,1)==completeParticleId(i));
end
nPore = height(poreInfo);
nParticle = height(particleInfo);
nThroat = height(throatConn);

DEM_Radius = DEMInfo(index_list(:,1) ,4);
DEM_Pos = DEMInfo(index_list(:,1) ,1:3);

formDrag = zeros(nParticle, 3);
skinFriction = zeros(nParticle, 3);
surfaceArea = zeros(nParticle, 2);
corrector = zeros(nParticle, 2);
%% calculate form drag
for i = id_complete
    id_interface =  find(spInterface.index_1==particleInfo.id(i) & spInterface.index_2<=nPore);
    id_pores = spInterface.index_2(id_interface);
    %
    Fpx = 0; Fpy = 0; Fpz = 0; 
    for j = 1:length(id_pores)
        porePressure = P_pore(id_pores(j));
        load(sprintf('%sthroat_solid_pore/%i__%i.mat', CasePath, id_interface(j)-1, particleInfo.id(i)), 'T')
        [Fpx_tmp, Fpy_tmp, Fpz_tmp, area_tmp] = integrateSurfacePressure(porePressure, particleInfo.xc(i), particleInfo.yc(i), particleInfo.zc(i), particleInfo.radius(i), T.x*dx, T.y*dx, T.z*dx, dx, []); 
        Fpx = Fpx+Fpx_tmp;Fpy = Fpy+Fpy_tmp;Fpz = Fpz+Fpz_tmp;
        surfaceArea(i, 1) = surfaceArea(i, 1) + area_tmp;
    end
    surfaceArea(i, 2) = 4*pi*(particleInfo.radius(i))^2;
    corrector(i) = surfaceArea(i, 2)./surfaceArea(i, 1);
    formDrag(i, 1:3) = [Fpx, Fpy, Fpz]*corrector(i);
    
end

%% calculate skin friction
for i = 1:nThroat
    % calculate the total viscous force
    id_pi = throatConn.poreI(i);
    id_pj = throatConn.poreJ(i);
    if id_pi>0&&id_pj>0
        vec_pij = ([poreInfo.xc(id_pj)-poreInfo.xc(id_pi), poreInfo.yc(id_pj)-poreInfo.yc(id_pi), poreInfo.zc(id_pj)-poreInfo.zc(id_pi)]);
        vec_pij = vec_pij/norm(vec_pij);
        Ff = (P_pore(id_pi)-P_pore(id_pj))*throatConn.radius(i)^2*pi*vec_pij;
        
        id_interface_i = find(spInterface.index_1==id_pi&spInterface.index_2>nPore);
        id_interface_j = find(spInterface.index_1==id_pj&spInterface.index_2>nPore);
        id_particles = intersect(spInterface.index_2(id_interface_i), spInterface.index_2(id_interface_j));
        
        A_particle = zeros(size(id_particles));
        for j = 1:length(id_particles)
            A_particle(j) = spInterface.Area(spInterface.index_1==id_pi&spInterface.index_2==id_particles(j)) + spInterface.Area(spInterface.index_1==id_pj&spInterface.index_2==id_particles(j));
        end
        A_total = sum(A_particle);
        
        for j = 1:length(id_particles)
            skinFriction(id_particles(j)-nPore, :) = skinFriction(id_particles(j)-nPore, :) + Ff*A_particle(j)/A_total;
        end
    end
end

formDrag = formDrag*[0,0,1;1,0,0;0,1,0];
skinFriction = skinFriction*[0,0,1;1,0,0;0,1,0];

Fd = formDrag + skinFriction;
Fd_abs = vecnorm(Fd, 2, 2);
formDrag_abs = vecnorm(formDrag, 2, 2);
skinFriction_abs = vecnorm(skinFriction, 2, 2);
sz = [length(Fd_abs), 6];
varNames = {'id_p', 'Rp', 'Fd', 'formDrag', 'skinFriction', 'stokesDrag'};
varTypes = {'double', 'double', 'double', 'double', 'double', 'double'};
PN_force = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
PN_force.id_p = index_list(:,1);
PN_force.Rp = particleInfo.radius;
PN_force.Fd = Fd;
PN_force.formDrag = formDrag;
PN_force.skinFriction = skinFriction;
PN_force.stokesDrag = 6*pi*particleInfo.radius;

