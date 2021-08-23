B=importdata('/home/rochi/work/monteCarloSolver/start_merge_MCDC/MCDC_Simulator_public/docs/scheme_files/PGSE_permeability_paper_scheme.scheme');
Bdata=B.data;
%%


A=importdata('perm_paper_cylinder_pr1_DWI.txt');
bval=zeros(length(Bdata),1);
count=1;
Gamma=3.3075*10^8;
for i=1:length(Bdata)
    bval(count)=(Gamma^2)*(Bdata(i,6)^2)*(Bdata(i,4)^2)*(Bdata(i,5)-Bdata(i,6)/3)*10^(-6);
    count=count+1;
end
Result=[Bdata bval A./10000];
ResultD1x=Result(1:10:500,:);
ResultD1y=Result(2:10:500,:);
ResultD2x=Result(3:10:500,:);
ResultD2y=Result(4:10:500,:);
ResultD3x=Result(5:10:500,:);
ResultD3y=Result(6:10:500,:);
ResultD4x=Result(7:10:500,:);
ResultD4y=Result(8:10:500,:);
ResultD5x=Result(9:10:500,:);
ResultD5y=Result(10:10:500,:);
figure(1)
%subplot(3,1,1)
semilogy(ResultD1y(:,4).^2,ResultD1y(:,9)./max(ResultD1y(:,9)),'.-','Color','k')
hold on;
semilogy(ResultD2y(:,4).^2,ResultD2y(:,9)./max(ResultD2y(:,9)),'.-','Color','k')
hold on;
semilogy(ResultD3y(:,4).^2,ResultD3y(:,9)./max(ResultD3y(:,9)),'.-','Color','k')
hold on;
semilogy(ResultD4y(:,4).^2,ResultD4y(:,9)./max(ResultD4y(:,9)),'.-','Color','k')
hold on;
semilogy(ResultD5y(:,4).^2,ResultD5y(:,9)./max(ResultD5y(:,9)),'.-','Color','k')
ylim([0.08 1.1])
xlim([0 1.1])
pbaspect([2 1 1])
xlabel('g^2 [T^2/m^2]') 
ylabel('Signal/Signal_o') 
title('p=0.0002')