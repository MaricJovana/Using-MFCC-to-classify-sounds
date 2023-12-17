clear all; close all; clc;
%radni folder mora da bude otvorena Velika baza da bi program radio
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
file_name=sprintf('broj_1_1.wav');               
        [y fs]=audioread(file_name);
        Tw=ceil(((length(y)/fs)*1000)/w);
        Ts=0.5*Tw;
        [MFCCs,FBEs,frames]=mfcc(y, fs, Tw, Ts, alpha, hamming, R, M, C(1), L ); %promeniti C ovde
        MFCC_vektor=reshape(MFCCs,1,size(MFCCs,1)*size(MFCCs,2));
        MFCC_svi_vektori(1,:)=MFCC_vektor;
%ucitavanje ostalih signala
for br1=1:N(1)
    for br2=1:N(2)
        file_name=sprintf('broj_%d_%d.wav', br1, br2);               
        [y fs]=audioread(file_name);
        Tw=ceil(((length(y)/fs)*1000)/w);
        Ts=0.5*Tw;
        [MFCCs,FBEs,frames]=mfcc(y, fs, Tw, Ts, alpha, hamming, R, M, C(1), L ); %promeniti C i ovde
        MFCC_vektor=reshape(MFCCs,1,size(MFCCs,1)*size(MFCCs,2));
        if(length(MFCC_vektor)>length(MFCC_svi_vektori)) %veci broj prozora nego fiksan
            MFCC_svi_vektori(((br1-1)*N(2)+br2),:)=MFCC_vektor(1:length(MFCC_svi_vektori));
        elseif(length(MFCC_vektor)<length(MFCC_svi_vektori)) %manji broj prozora nego fiksan
            ze =zeros(1,length(MFCC_svi_vektori(1,:))); 
            ze(1:length(MFCC_vektor))=MFCC_vektor(1:length(MFCC_vektor));
            MFCC_svi_vektori(((br1-1)*N(2)+br2),:)=ze;
        else %isti broj prozora kao fiksan
            MFCC_svi_vektori(((br1-1)*N(2)+br2),:)=MFCC_vektor;
        end
    end
end

%EUKLIDSKA DISTANCA
for br3=1:N(3)
    for br4=1:N(3)
        D(br3,br4)=norm(MFCC_svi_vektori(br4,:)-MFCC_svi_vektori(br3,:));
    end
end
%figure, surf(D); colorbar; view(90,90); title('Euklidska distanca','EdgeColor','none');colormap jet;
figure, imagesc(D);axis('xy');colorbar;title('Euklidska distanca'); colormap jet;

%ODLUCIVANJE
prag = 80; %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~menjati prag ovde
%alternativa sa for petljama
%{
for br3=1:N(3)
    for br4=1:N(3)
        if(D(br3,br4)<prag) 
            D_odl(br3,br4)=1;
        else
            D_odl(br3,br4)=0;
        end
    end
end
%}
D_odl=D<prag;
%figure, surf(D); colorbar; view(90,90); title('Prepoznavanje');
figure, imagesc(D_odl); colorbar; view(90,90); title('Prepoznavanje');colormap jet;
%idealna matrica prepoznavanja
D_ideal = zeros(N(3),N(3));
D_ideal(1:N(2),1:N(2))=1;
for br5=1:N(1)-1
    D_ideal((N(2))*br5+1:((br5+1)*N(2)),N(2)*br5+1:((br5+1)*N(2)))=1;
end
%figure, surf(D_ideal); colorbar; view(90,90); title('Idealan rezultat');
figure, imagesc(D_ideal); colorbar; view(90,90); title('Idealan rezultat');
%procenat
for br3=1:N(3)
    for br4=1:N(3)
        D_procenat(br3,br4)=abs(D_odl(br3,br4)-D_ideal(br3,br4));
    end
end
br_gresaka=sum(sum(D_procenat));
procenat_greske = (br_gresaka/(N(3))^2)*100; %greske
procenat_uspesnosti = 100-procenat_greske;
display(['Procenat tacnih odluka ', num2str(procenat_uspesnosti) '%']);
sound(y,fs); %samo cujem kad je gotov program

%zapisane vrednosti posle vise pokretanja programa
p=[70 75 80 85 90 100];
pr=[86.1899 87.3515 87.8817 87.3959 85.6857 78.1216];
figure,subplot(3,1,1), plot(p,pr); title('za 3 MFCC koeficijenata'); 
p=[100 110 120 130 140];
pr=[84.8473 87.3886 88.4126 85.7643 77.9558];
subplot(3,1,2), plot(p,pr); title('za 6 MFCC koeficijenata'); 
p=[100 110 120 130 150 160 170 180];
pr=[80.2089 80.2902 80.4796 80.9847 83.5298 85.0453 85.6619 83.6547];
subplot(3,1,3), plot(p,pr); title('za 12 MFCC koeficijenata'); 
sgtitle('Zavisnost procenta uspesnosti prepoznavanja reci od vrednosti praga');
