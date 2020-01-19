GetDatainMediumUpdatedConcUpdate
load('ParticleTypeKey.mat')
chemistries={'ZnO','Ag', 'CeO2', 'TiO2' ,'CuO','Fe2O3','Fe3O4'};
count=1;

for i=1:length(Data)
    %Checks if complete size data was provided
    if isnan(Data(i,4)) || Data(i,4)==-9999 ||Data(i,5)==9999
        continue
    end
    %checks if complete concentration data was provided
    if isnan(Data(i,9))
        continue
    end 
%checks to be sure proper chemistry is used
    if sum(strcmp(chemistries,ParticleTypeKey{Data(i,1)}))==0
        continue
    end
    
    if Data(i,12)==1
        ConcentrationData(count,1)=i;
        ConcentrationData(count,2)=Data(i,8);
        ConcentrationData(count,3)=(Data(i,9))/Data(i,8);
        ConcentrationData(count,4)=particlecalc(Data(i,8),Data(i,3),Data(i,1));
        ConcentrationData(count,5)=(ConcentrationData(count,3)^2+...
            (3*(Data(i,3)-Data(i,4))/Data(i,3))^2)^(1/2);
        ConcentrationData(count,6)=Surfacalc(Data(i,8),Data(i,3),Data(i,1));
        ConcentrationData(count,7)=(ConcentrationData(count,3)^2+...
            ((Data(i,3)-Data(i,4))/Data(i,4))^2)^(1/2);    
        count=count+1;
    end
end
clear count i
ConcentrationData(ConcentrationData(:,4)==-9999,:)=[];
ConcentrationData(ConcentrationData(:,3)==0,:)=[];
count=1;
for i=1:length(Data(ConcentrationData(:,1),7))
    %Captures if there is a number in the tox description
    [begin, stopper]=regexp(ToxTypeKey{Data(ConcentrationData(i,1),7)},...
        '\d.\.?\d*');
    %Looks to see if there is only one number in the key and then if the
    %number is 50. If so the data is added to the EP50 data
    if length(begin)==1 && str2double(...
            ToxTypeKey{Data(ConcentrationData(i,1),7)}(begin:stopper))==50
        EP50data(count,:)=ConcentrationData(i,:);
        count=count+1;
    end
    %checking if 50 is in the interval of toxicity 
    if length(begin)==2 && str2double(ToxTypeKey{...
            Data(ConcentrationData(i,1),7)}(begin(1):stopper(1)))+...
            str2double(ToxTypeKey{Data(ConcentrationData(i,1),7)}...
            (begin(2):stopper(2)))>=50 &&str2double(ToxTypeKey{...
            Data(ConcentrationData(i,1),7)}(begin(1):stopper(1)))-...
            str2double(ToxTypeKey{Data(ConcentrationData(i,1),7)}...
            (begin(2):stopper(2)))<=50
        EP50data(count,:)=ConcentrationData(i,:);
        count=count+1;
    end
end
count=1;
%Selects Records that have mortality as their endpoint
Regexphrases=["mortality","immo","lum"];
for i=1:length(EP50data)
    for j=1:length(Regexphrases)
        if ~isempty(regexp(CellDat{Data(EP50data(i,1),13), 22}, Regexphrases(j),'once'))
            Outputhold(count,:)=EP50data(i,:);
            count=count+1;
        end 
    end 
end
clear EP50data
EP50data=Outputhold;
clear Outputhold
coatingregex=["#N/A", "(coated)$", "cit","prote", "PVP","doped",...
    "alka","BPEI"];

for i=1:length(EP50data)
    for j=1:length(coatingregex)
        if ~isempty(regexp(CellDat{Data(EP50data(i,1),13),5}, coatingregex(j)))
            EP50data(i,8)=j;
        end
    end
end

Shapeaccount1=EP50data;
Shapeaccount2=EP50data;
Shapeaccount1(ismember(Data(Shapeaccount1(:,1),2),...
    [1,17,18,21,35,40,41]),:)=[];
Shapeaccount2(ismember(Data(Shapeaccount2(:,1),2),...
    [21,40]),:)=[];
Shapeaccount1(Data(Shapeaccount1(:,1),10)==1,:)=[];
Shapeaccount2(Data(Shapeaccount2(:,1),10)==1,:)=[];

csvwrite('RData\ITMCompShape.csv',Shapeaccount1)
csvwrite('RData\ITMComp.csv',Shapeaccount2)

clear ConcentrationData EP50data i begin stopper count %Shapeaccount1



%Allowing Data with no size error
count=1;
for i=1:length(Data)
    %checks to make sure there is a size value given and not one of the
    %indistinct error points
    if isnan(Data(i,3)) || Data(i,4)==-9999 ||Data(i,5)==9999
        continue
    end
    if isnan(Data(i,9))
        continue
    end 
%checks to be sure proper chemistry is used
    if sum(strcmp(chemistries,ParticleTypeKey{Data(i,1)}))==0
        continue
    end
    
    if Data(i,12)==1
        ConcentrationData(count,1)=i;
        ConcentrationData(count,2)=Data(i,8);
        ConcentrationData(count,3)=(Data(i,9))/Data(i,8);
        ConcentrationData(count,4)=particlecalc(Data(i,8),Data(i,3),Data(i,1));
        ConcentrationData(count,5)=(ConcentrationData(count,3));
        ConcentrationData(count,6)=Surfacalc(Data(i,8),Data(i,3),Data(i,1));
        ConcentrationData(count,7)=(ConcentrationData(count,3));
        count=count+1;
    end
end
clear count i
ConcentrationData(ConcentrationData(:,4)==-9999,:)=[];
ConcentrationData(ConcentrationData(:,3)==0,:)=[];
count=1;
for i=1:length(Data(ConcentrationData(:,1),7))
    [begin, stopper]=regexp(ToxTypeKey{Data(ConcentrationData(i,1),7)},...
        '\d.\.?\d*');
    if length(begin)==1 && str2double(...
            ToxTypeKey{Data(ConcentrationData(i,1),7)}(begin:stopper))==50
        EP50data(count,:)=ConcentrationData(i,:);
        count=count+1;
    end
    %checking if 50 is in the interval of toxicity 
    if length(begin)==2 && str2double(ToxTypeKey{...
            Data(ConcentrationData(i,1),7)}(begin(1):stopper(1)))+...
            str2double(ToxTypeKey{Data(ConcentrationData(i,1),7)}...
            (begin(2):stopper(2)))>=50 &&str2double(ToxTypeKey{...
            Data(ConcentrationData(i,1),7)}(begin(1):stopper(1)))-...
            str2double(ToxTypeKey{Data(ConcentrationData(i,1),7)}...
            (begin(2):stopper(2)))<=50
        EP50data(count,:)=ConcentrationData(i,:);
        count=count+1;
    end
end
count=1;
for i=1:length(EP50data)
    for j=1:length(Regexphrases)
        if ~isempty(regexp(CellDat{Data(EP50data(i,1),13), 22}, Regexphrases(j),'once'))
            Outputhold(count,:)=EP50data(i,:);
            count=count+1;
        end 
    end 
end
clear EP50data
EP50data=Outputhold;
clear Outputhold
for i=1:length(EP50data)
    for j=1:length(coatingregex)
        if ~isempty(regexp(CellDat{Data(EP50data(i,1),13),5}, coatingregex(j)))
            EP50data(i,8)=j;
        end
    end
end
Shapeaccount3=EP50data;
Shapeaccount3(ismember(Data(Shapeaccount3(:,1),2),...
    [21,40]),:)=[];
Shapeaccount3(Data(Shapeaccount3(:,1),10)==1,:)=[];
csvwrite('RData\ITMIncomp.csv',Shapeaccount3)

