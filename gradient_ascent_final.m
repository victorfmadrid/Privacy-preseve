%% variables
clear all
clc
ps = [1/5 2/5 2/5];
pz = [1/4 3/4];
len_y = 9;
len_u = 11;
len_z = 2;
psz = [];%ps x pz
%% Find P(s,z)
for i = 1:length(ps)
    for j = 1:length(pz)
        psz(i,j) = ps(i) * pz(j);
    end
end
%% Find P(y|s,z=1)
len_s = 3;
py_sz1= []; %probabilities of p(s|y(i))
dist_1 = [2 7 5];
x = 0:len_y-1;
for i = 1:length(ps)
    reg = poisspdf(x,dist_1(i));
    reg = reg / sum(reg);
    py_sz1(i,:) = reg;%each row is a different lambda value for poisson dist
end
py_sz1 = py_sz1';%val of y x val of s

%% Find P(y|s,z=0)
py_sz0= []; %probabilities of p(s|y(i))
dist_1 = [5 4 8];
x = 0:len_y-1;
for i = 1:length(ps)
    reg = poisspdf(x,dist_1(i));
    reg = reg / sum(reg);
    py_sz0(i,:)= reg;%each row is a different lambda value for poisson dist
end
py_sz0 = py_sz0';%val of y x val of 
py_sz = cell(1,2);
py_sz{1,1} = py_sz0;
py_sz{1,2} = py_sz1;
%% Find P(y,s,z)
psz = psz';
pysz0 = [];
for i=1:len_y
    for j = 1:length(ps)
            pysz0(i,j) = py_sz0(i,j) * psz(1,j);
    end
end

pysz1 = [];
for i=1:len_y
    for j = 1:length(ps)
            pysz1(i,j) = py_sz0(i,j) * psz(2,j);
    end
end
pysz = cell(1,2);
pysz{1,1} = pysz0;
pysz{1,2} = pysz1;
%pysz is a ZxYxS
%% P(Y)
py = sum(pysz0,2) + sum(pysz1,2);%gives you p(Y)
%% Intitial mapping of P(u|y)
pu_y = []; %probabilities of p(s|y(i))
dist_m = [3,7,4,8,7,2,8,9,5,8,6];
x = 0:len_u-1;
for i = 1:1:len_y
    reg = poisspdf(x,dist_m(i));
    reg = reg / sum(reg);
    pu_y(:,i) = reg;
end
%this is done columnwise sum should be 1
%% Intial mapping of pu
pu = [];
for i = 1:len_u
    sum_box = [];
    for j = 1:len_y
        temp = pu_y(i,j) * py(j);
        sum_box = [sum_box,temp];
    end
    pu(i) = sum(sum_box);
end
%% Mappings of P(S,Y) and P(Z,Y)
pzy = [];
pzy(:,1) = sum(pysz0,2);
pzy(:,2) = sum(pysz1,2);

psy = pysz0 + pysz1;

%% Intial mapping of p(z|u)
pz_u =[];
pzy = pzy';
for i = 1:length(pz)
    sum_box2 = [];
    for j = 1:length(pu)
        sum_box = [];
        for k = 1:length(py)
            temp = pu_y(j,k)*pzy(i,k);
            sum_box = [sum_box, temp];
        end 
        temp2 = sum(sum_box) / pu(j);
        sum_box2 = [sum_box2,temp2];
    end
     pz_u(i,:)= sum_box2;
end
%% Initial mappig of p(s|u)
ps_u = []; 
psy = psy';
for i = 1:len_s
    sum_box2 = [];
    for j = 1:length(pu)
        sum_box = [];
        for k = 1:length(py)
            temp = pu_y(j,k)*psy(i,k);
            sum_box = [sum_box, temp];
        end 
        temp2 = sum(sum_box) / pu(j);
        sum_box2 = [sum_box2,temp2];
    end
     ps_u(i,:)= sum_box2;
end

%% Intial mapping for LAMBDA
lmda = [3.17,8.14,7.89,8.52,5.05,6.35,9.50,4.43,0.60,8.66,6.31];
%% Variables
alpha = 2;
beta = 7;
rho = 1;
%% P(S|Yi)
ps_y = [];
for i=1:len_s
    for j=1:len_y
        ps_y(i,j) = psy(i,j)/py(j);
    end
end
%% P(Z|Yi)
pz_y = [];
for i=1:length(pz)
    for j=1:len_y
        pz_y(i,j) = pzy(i,j)/py(j);
    end
end
%% Before setting up loop _________________________________________________
initial_puy = pu_y;
intial_pu = pu;
step = .01;
quick_graph = [];
for B=1:40
    %% Finding I(S;U)
    isu =[];
    for i=1:len_u
        row = [];
        for j=1:len_y
            temp1 = [];
            for k=1:len_s
                temp2 = ps_y(k,j) * (log( ps_y(k,j) /ps_u(k,i)));
                temp1 = [temp1,temp2];
            end
            entry = -1* sum(temp1) *py(j);
            row = [row,entry];
        end
        isu(i,:) = row;
    end
    %% Finding Alpha*I(Z;U)
    izu =[];
    for i=1:len_u
        row = [];
        for j=1:len_y
            temp1 = [];
            for k=1:len_z
                temp2 = pz_y(k,j) * (log( pz_y(k,j) /pz_u(k,i)));
                temp1 = [temp1,temp2];
            end
            entry = alpha*sum(temp1) *py(j);
            row = [row,entry];
        end
        izu(i,:) = row;
    end
    %% Finding lambda(u)*delta(u)
    delta1 = [];
    for i=1:length(lmda)
        for j = 1:len_y
            delta1(i,j) = -1*lmda(i) *py(j);
        end
    end

    %% Finding rho/2*delta(u)^2
    delta2 = [];
    for i=1:len_u
        for j = 1:len_y
            delta2(i,j) = rho * py(j) * (pu(i) - pu_y(i,j)*py(j));
        end
    end
    %% Putting it together to find P(U|Y)^t

    puy_t = [];
    for i=1:len_u
        for j = 1:len_y
            puy_t(i,j) = isu(i,j) +izu(i,j) + delta1(i,j) + delta2(i,j);
        end

    end
    
    Puy_t1 =[];
    for i =1:len_u
        for j= 1:len_y
            Puy_t1(i,j) = initial_puy(i,j) - step*puy_t(i,j);
        end
    end
    initial_puy = Puy_t1;
    quick_graph = [quick_graph,Puy_t1(1,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %this is the part where we find p(u) using P(U|Y)^t
    %% Finding Finding rho/2*delta(u)^2

    delta3 = [];
    for i = 1:len_u
        sum_box = [];
        for j = 1:len_y
            temp = Puy_t1(i,j) * py(j);
            sum_box = [sum_box,temp];
        end
        delta3(i) = -1* rho *( pu(i) - sum(sum_box));
    end
    %% putting it together to find p(u)t

    pu_t=[];
    for i = 1:len_u
        pu_t(i) = lmda(i) + delta3(i);
    end
    Pu_t1 =[];
    for i = 1:len_u
        Pu_t1 = intial_pu - step* pu_t(i);
    end
    intial_pu = Pu_t1;


    
    %% Finding new lmba
    new_lambda =[];
    for i=1:len_u
        add = [];
        for j =1:len_y
            temp = Puy_t1(i,j) * py(j);
            add = [add,temp];
        end
        new_lambda(i) = lmda(i) + rho*(Pu_t1(i) - sum(add));
    end

    %% update P(z|u) 
    pz_u =[];
    for i = 1:len_z
        sum_box2 = [];
        for j = 1:len_u
            sum_box = [];
            for k = 1:len_y
                temp = Puy_t1(j,k)*pzy(i,k);
                sum_box = [sum_box, temp];
            end 
            temp2 = sum(sum_box) / Pu_t1(j);
            sum_box2 = [sum_box2,temp2];
        end
         pz_u(i,:)= sum_box2;
    end
    %% updata p(s|u)
    ps_u = []; 
    for i = 1:len_s
        sum_box2 = [];
        for j = 1:len_u
            sum_box = [];
            for k = 1:len_y
                temp = Puy_t1(j,k)*psy(i,k);
                sum_box = [sum_box, temp];
            end 
            temp2 = sum(sum_box) / Pu_t1(j);
            sum_box2 = [sum_box2,temp2];
        end
         ps_u(i,:)= sum_box2;
    end
end
%% Plot curve
figure()
plot(1:B,abs(quick_graph));


    



