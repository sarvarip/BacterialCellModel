function F = simulation_Algar(t, Y, L, M, alfaplus, alfaminus, beta)

% L = 50;
% M = 3;
% alfaplus = 10^-4;
% alfaminus = 200;
% beta = ones(L+1,1); %everything shifted by one, since there is no idx zero in Matlab

%Settings - to compare

% L = 1000;
% M = 1;
% Tini = 1;
% alfaplus = 10^-6/Tini; %Vini = Kini*M*R
% alfaminus = 0;
% R = 10^6; %a lot, and M is small -> neglect consumption of ribosomes
% Tel = 4; %100000 edge case
% beta = 1/Tel*ones(L+1,1); %everything shifted by one, since there is no idx zero in Matlab

%G is Y(1), Y(2) is Y0, Y(3) is Y1 .. so everything shifted by 2
%Beta0 is beta(1), Beta1 is beta(2).. so everything shifted by 1

F = zeros(L+2,1);

F(1) = -M*alfaplus*Y(1)*(1-Y(2)/M)+alfaminus*Y(2)+beta(L+1)*Y(L+2); 
F(2) = M*alfaplus*Y(1)*(1-Y(2)/M)-alfaminus*Y(2)-beta(1)*Y(2)*(1-Y(3)/M);

for i=3:L+1
    F(i) = beta(i-2)*Y(i-1)*(1-Y(i)/M)-beta(i-1)*Y(i)*(1-Y(i+1)/M);
end

F(L+2) = beta(L)*Y(L+1)*(1-Y(L+2)/M)-beta(L+1)*Y(L+2);