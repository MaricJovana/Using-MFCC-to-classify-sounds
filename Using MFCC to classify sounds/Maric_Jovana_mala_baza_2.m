clear all; close all; clc;

alpha=0.97;      
R=[300 3700 ];  
M=30;            
C=[3,6,12];            
L=22;            
hamming=@(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));
N=[5,205,1025];
w=15;

%UCITAVANJE SIGNALA, RACUNANJE MFCC
%prvi snimak radim odvojeno da bi se odredio fiksan broj prozora za sve
file_name=sprintf('broj1_1.wav');               
        [y fs]=audioread(file_name);
        Tw=ceil(((length(y)/fs)*1000)/w);
        Ts=0.5*Tw;
        [MFCCs,FBEs,frames]=mfcc(y, fs, Tw, Ts, alpha, hamming, R, M, C(1), L ); %promeniti C ovde
        MFCC_vektor=reshape(MFCCs,1,size(MFCCs,1)*size(MFCCs,2));
        MFCC_svi_vektori(1,:)=MFCC_vektor;
%ucitavanje ostalih signala
for br1=1:N(1)
    for br2=1:N(1)
        file_name=sprintf('broj%d_%d.wav', br1, br2);               
        [y fs]=audioread(file_name);
        Tw=ceil(((length(y)/fs)*1000)/w);
        Ts=0.5*Tw;
        [MFCCs,FBEs,frames]=mfcc(y, fs, Tw, Ts, alpha, hamming, R, M, C(1), L ); %promeniti C i ovde
        MFCC_vektor=reshape(MFCCs,1,size(MFCCs,1)*size(MFCCs,2));
        if(length(MFCC_vektor)>length(MFCC_svi_vektori)) %veci broj prozora nego fiksan
            MFCC_svi_vektori(((br1-1)*N(1)+br2),:)=MFCC_vektor(1:length(MFCC_svi_vektori));
        elseif(length(MFCC_vektor)<length(MFCC_svi_vektori)) %manji broj prozora nego fiksan
            ze =zeros(1,length(MFCC_svi_vektori(1,:))); 
            ze(1:length(MFCC_vektor))=MFCC_vektor(1:length(MFCC_vektor));
            MFCC_svi_vektori(((br1-1)*N(1)+br2),:)=ze;
        else %isti broj prozora kao fiksan
            MFCC_svi_vektori(((br1-1)*N(1)+br2),:)=MFCC_vektor;
        end
    end
end

%EUKLIDSKA DISTANCA
for br3=1:N(1)^2
    for br4=1:N(1)^2
        D(br3,br4)=norm(MFCC_svi_vektori(br4,:)-MFCC_svi_vektori(br3,:));
    end
end
%figure, surf(D); colorbar; view(90,90); title('Euklidska distanca');
figure, imagesc(D);axis('xy');colorbar;title('Euklidska distanca'); colormap jet;

%ODLUCIVANJE
prag = 76; %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~menjati prag ovde
% 100 za 6 koeficijenata, 120 za 12 koeficijenata, 76 za 3 koeficijenta
for br3=1:N(1)^2
    for br4=1:N(1)^2
        if(D(br3,br4)<prag) 
            D_odl(br3,br4)=1;
        else
            D_odl(br3,br4)=0;
        end
    end
end

%figure, surf(D); colorbar; view(90,90); title('Prepoznavanje');
figure, imagesc(D_odl); colorbar; view(90,90); title('Prepoznavanje');
%idealna matrica prepoznavanja
D_ideal = zeros(N(1)^2,N(1)^2);
D_ideal(1:N(1),1:N(1))=1;
for br5=1:N(1)-1
    D_ideal((N(1))*br5+1:((br5+1)*N(1)),N(1)*br5+1:((br5+1)*N(1)))=1;
end
%figure, surf(D_ideal); colorbar; view(90,90); title('Idealan rezultat');
figure, imagesc(D_ideal); colorbar; view(90,90); title('Idealan rezultat');
%procenat
for br3=1:N(1)^2
    for br4=1:N(1)^2
        D_procenat(br3,br4)=abs(D_odl(br3,br4)-D_ideal(br3,br4));
    end
end
br_gresaka=sum(sum(D_procenat));
procenat_greske = (br_gresaka/(N(1))^4)*100; %greske
procenat_uspesnosti = 100-procenat_greske;
display(['Procenat tacnih odluka ', num2str(procenat_uspesnosti) '%']);
sound(y,fs); %samo cujem kad je gotov program