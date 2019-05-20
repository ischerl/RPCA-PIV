%Function to apply vector interpolation to specified data set
%Owen Williams

%Add possible filter to exclude the interpolation of vectors missing a
%certain number of neighbours

function InterpVectors

close all 
clear all

W = NaN;
Sn = NaN;
Info = NaN;

[File,Folder]=uigetfile('*.mat', 'Pick data to calculate missing vectors: ');
load(strcat(Folder, File));

numData = size(U);
MissingVecs = zeros(numData);
MissingVecs(isnan(U)) = 1;
numData = numData(3);
PercentMissing = zeros(numData,1);
tags = NaN;

while 1
    fprintf(1,'\n')
    disp('Select method to interpolate missing vectors:  ')
    disp('  0:   Exit (without saving)                   ')
    disp('  1:   Linear Interpolation                    ')
    disp('  2:   Weighted Interpolation                  ')
    disp('  99:  Save Data                               ')
    
    fprintf(1,'\n')
    in = input('Input your option here: ');
    
    switch in
        case 0
            clc
            return
            
        case 1
            if isnan(tags)
                for i = 1:numData
                   PercentMissing(i) = sum(sum(MissingVecs(:,:,i)))/numel(MissingVecs(:,:,i))*100; 
                end
                
%                 gui_active(1); % will add an abort button
%                 h = progressbar( [],0,'Calculating Missing Vectors (Linear)' );
                
                for i = 1:numData
                    [U(:,:,i),V(:,:,i)]=naninterp(U(:,:,i),V(:,:,i),'linear');
%                     h = progressbar( h,1/numData );
%                     if ~gui_active
%                         progressbar( h,-1 );
%                         return;
%                     end
                end
                
%                 progressbar( h,-1 );
                tags = '-IntLin';
            else
                disp('Vectors are already interpolated!');
            end
            
         case 2
            if isnan(tags)
                for i = 1:numData
                   PercentMissing(i) = sum(sum(MissingVecs(:,:,i)))/numel(MissingVecs(:,:,i))*100; 
                end
                
%                 gui_active(1); % will add an abort button
%                 h = progressbar( [],0,'Calculating Missing Vectors (Weighted)' );
                
                for i = 1:numData
                    [U(:,:,i),V(:,:,i)]=naninterp(U(:,:,i),V(:,:,i),'weighted');
%                     h = progressbar( h,1/numData );
%                     if ~gui_active
%                         progressbar( h,-1 );
%                         return;
%                     end
                end
                
%                 progressbar( h,-1 );
                tags = '-IntWei';
            else
                disp('Vectors are already interpolated!');
            end

        case 99
            if sum(PercentMissing)==0
                disp('No vectors have been interpolated yet!')
            else
                %newFile = input('Save as (default: Original Filename-tags): ');
%                 if isempty(newFile)
%                     File = strcat(File(1:end-4), tags);
%                 else
%                     File=newFile;
%                 end
                File = 'vel_int';
                
%                 if isnan(W)
%                     if isnan(Sn)
%                         save(strcat(Folder,File), 'U','V', 'source', 'PercentMissing')
%                     else
%                         save(strcat(Folder,File), 'X','Y','U','V', 'Sn','Info', 'source', 'PercentMissing')
%                     end
%                 else
%                     save(strcat(Folder,File), 'X','Y','U','V','W', 'source', 'PercentMissing')
%                 end
                save(strcat(Folder,File), 'U','V','PercentMissing')
                
                disp(['File saved as : ' File])
            
            end
 
        otherwise
            fprintf(1, 'Invalid input\n');
            
        input('\nHit Enter to Proceed')
        
    end
end

end