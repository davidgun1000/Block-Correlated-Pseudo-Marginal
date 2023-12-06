%adjust data
load('DiffusionDataDailyRaw.mat');
Raw = [Aero,Agric,Autos,Banks,Beer,BldMt,Books,Boxes,BusSv,Chems,Chips,Clths,Cnstr,Coal,Drugs,ElcEq,FabPr,Fin,Food,Fun,Gold,Guns,Hardw,Hlth,Hshld,Insur,LabEq,Mach,Meals,MedEq,Mines,Oil,Other,Paper,PerSv,RlEst,Rtail,Rubbr,Ships,Smoke,Soda,Softw,Steel,Telcm,Toys,Trans,Txtls,Util,Whlsl];
Raw = Raw(1:100,:);
y=Raw' - mean(Raw',2);
save('DiffusionDataDaily.mat','y');
