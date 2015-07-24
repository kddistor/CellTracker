%CELLTRACER_TEMP
%   Script for processing cell image data.
%   THIS SCRIPT IS TO BE EDITTED FOR EACH PROCESSING RUN
%
%   Best practice recommentdations.
%   1.  Ensure that your MATLAB workding directory is on your local machine
%       (or connect high speed, e.g. USB 3.0).  It will need to be able to
%       accommodate tens (10) of gigabytes.
%   2.  Ensure that you set the temporary save directory to be on either
%       your local drive or a high speed (e.g. USB 3.0) connection.
%       WARNING, the temporary saves may accumulate up to hundreds (100) of
%       gigabytes and many local drives will be insufficient.
%   3.  If setting this prodecure up for the first time on your machine,
%       first run the provided block of test code to evaluate the function
%       of the processing code (uncomment toward end of this script).
%   
%   Example usage procedure:
%   1.  Copy data to be processed to local drive (i.e. the hard drive of
%       your computer or an external drive connected with USB 3.0).  The
%       procedure will operate much faster if the data is local and any
%       copying from network locations has been performed independently.
%   2.  Specify 'Operation definitions' for your processing run.
%   3.  Specify 'Data definitions' for your processing run (name of source
%       file and a (new) save file name, which XY positions to process).
%   4.  Specify the background regions or values for each XY position.
%       a.  Open ND2 file in Elements Viewer (full version not required)
%       b.  Using the cursor in Elements, hover over the Upper Left and
%           Lower Right corner of a suitable background region, read the
%           coordinates, and type them into the 'bk(n).reg' definitions    
%           (as [ULx, ULy; LRx, LRy]). If the region is good for all time
%           points, set bk(n).dyn = true, otherwise, only consider the
%           first frame and set bk(n).dyn = false.
%           IF no suitable background region is available, provide a
%           conservative estimate of background values (e.g. via other
%           wells) as [Channel 1, Channel 2, Channel 3, etc.], and set
%           bk(n).fix = true.
%   5.  IMPORTANT:  Ensure that the Backup Metadata matches your
%       file.  Do this by right clicking on the image in Elements and
%       selecting 'Image Properties'.  The 'Image Fields' tab will include
%       data on the objective, pixel size (um/px), and exposure and camera
%       settings for each Channel.  Note that Fluorophores are never in the
%       ND2 Metadata and must be updated for each run.
%       a.  IF desired, backup values may be forced (regardless of values
%           in the ND2 metadata) by setting the bkenforce variable.
%


% --- Operation definitions --------------------------------------------
%Segmentation directions
op.chan = 2;                            %Channel to use for segmenting
op.local = 'nuc';                       %Indicate if nuclear/cytoplasmic
op.cname = {'CFP2','YFP2', 'RFP'};       %Channel names
%Nucleus identification definitions
%Max/Min Diameter               MinFormFactor     Performing diagnostic?
op.maxD = 40; op.minD = 12;     op.minF = 0.8;    op.diag = 0;
                %^ Modified with * sqrt(0.3) from previous script
op.storemask = false;    %Store cells masks?  (may accumulate memory)

%Flag to display (sparsely) samples from processing
sampdisp = false;
                
%Number of parallel workers desired
nW = 10;

%Define existing ND2 Reader, if present
r = [];  ri = [];   %No previous Reader


% --- Data definitions -------------------------------------------------
% %Source file name (full path)
% fname = ['L:\albeck\imageData\',...
%     '2015-06-04-TrkA-FIRE.nd2']; 
fname = ['K:\',...
    '2015-05-21-ERKTR-FIRE.nd2']; 

     
%Save file name
sname = '2015-05-21-ERKTR-FIRE';
% sname = 'K:\2015-05-19-5e-EKARFoxo-exosomes';
% 
% %Directory for temporary save data
% tname = 'E:\MatlabData';  %Obsolete. Only need for iman_celltrace_main

%XY Positions desired
xypos=25:32;

%Sampling time, in minutes
tsamp = 5;      

%Baseline value from camera (Andor Zyla 5.5 is typically 100)
bval = 100;     

%Channels expected to cross-talk (I: independent, D: dependent)
% xtalk.I = 'YFP'; xtalk.D = 'RFP';   %YFP is expected to bleed into RFP
xtalk.I = [];    xtalk.D = [];    %<- Use if not checking cross-talk

%Known XY shift during imaging (e.g. due to spike-in perturbation)
%   List multiple shifts relative to previous frame (not start of imaging)
% xyshift = [];   %No shift to report
xyshift.frame = [69, 76];   %First shifted frame
xyshift.dx = [19, -8];      %Magnitude of x shifts (x_after - x_before)
xyshift.dy = [-16, 2];      %Magnitude of y shifts (y_after - y_before)

%Define image background intensity regions
%   Open file in Elements Viewer and select background region
%   Define corners of a box (upper left; lower right), one per xy, and
%   designate each region as dynamic (good for every time point) or not.
%   Remember that indices cannot be zero (0), so upper left corner is 1,1.
clear bk
% s = 1;
% bk(s).reg =  [220, 17; 306, 120];       bk(s).fix = true;  s = s + 1;
% bk(s).reg =  [774, 710; 844, 785];      bk(s).fix = true;  s = s + 1;
% [bk(1:end).dyn] = deal(false);


% Note: for Fixed value background:
bk(1).reg = [561,136, 265]; bk(1).dyn = false; bk(1).fix = true;
bk(2:115) = bk(1);
%            ^Individual channel values in order (do not subtract baseline)


%% Backup MetaData for Imaging Experiment
%   IMPORTANT:  Verify accuracy of backup data before running
%               Comment unknown entries to prevent corruption
bkmd.obj.Desc       =   'Apo_20x';      %Description of Objective
bkmd.cam.Desc       =   'zyla5';      %Description of Objective
bkmd.obj.Mag        =   20;             %Magnification
bkmd.obj.WkDist     =   1;              %Working Distance (mm)
bkmd.obj.RefIndex   =   1;              %Index of Refraction
bkmd.obj.NA         =   0.75;           %Numerical Aperture

bkmd.cam.PixSizeX   = 	0.65;           %Pixel horizontal size (um)
bkmd.cam.PixSizeY   =   0.65;           %Pixel vertical size (um)
bkmd.cam.PixNumX    =   1280;           %Number of Pixels in X
bkmd.cam.PixNumY    =   1080;           %Number of Pixels in Y
bkmd.cam.BinSizeX   =   2;              %Number of Pixels per bin, X
bkmd.cam.BinSizeY   =   2;              %Number of Pixels per bin, Y

bkmd.exp.Channel    =   {'CFP2', 'YFP2', 'RFP2'};  %Name(s) of Channel(s)
bkmd.exp.Exposure   =   [200, 200, 400];        %Exposure times (ms)
bkmd.exp.ExVolt     =   [15, 16, 40];           %Relative voltages used
bkmd.exp.Filter     =   {'Filter_Turq', 'Filter_Venus', 'Filter_Cherry2'};
bkmd.exp.FPhore     =   {'CFP2', 'YFP2', 'RFP2'};
%   See iman_name for a list of valid Filter and Fluorophore Names

%IF desired, indicate backup fields to enforce (disregard ND2 file data)
%   In a Nx2 cell array, specify the category (obj, cam, exp) and field.
bkenforce = {'exp', 'Channel'};
% bkenforce = [];       %<- Use if not enforcing any fields


%%  TESTING SECTION
%  ---  FOR TESTING USE  ---  UNCOMMENT TO TEST PROCEDURE  ---
% fname = ['L:\albeck\Code\Image Processing\',...
%     'test1_crop01.nd2'];
% xypos = 1;    clear bk;
% xyshift = [];
% % xyshift.frame = [3, 5];  xyshift.dx = [-10; 20];  xyshift.dy = [-10; 20];
% sname = 'test1_crop01';        %tname = 'C:\Temp';
% nW = 1;     %Set nW to 1 to stop at breakpoints, higher to test parallel
% op.chan = 2;  op.local = 'cyt';  op.diag = false;  op.storemask = true;
% 
% s = 1;
% bk(s).reg =  [50, 1; 230, 45];          bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [1260, 1; 1280, 130];      bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [1000, 175, 275];          bk(s).fix = true;  s = s + 1;
% bk(s).reg =  [1, 1; 50, 35];            bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [930, 700; 1070, 830];     bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [10, 680; 105, 920];       bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [1130, 1; 1280, 45];       bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [1183, 992; 1280, 1080];	bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [1245, 1; 1280, 75];       bk(s).fix = false;  s = s + 1;
% bk(s).reg =  [1, 960; 80, 1080];        bk(s).fix = false;  s = s + 1;
% [bk(1:end).dyn] = deal(true);           
% bk(2).dyn = false;  bk(3).dyn = false;
% bk = bk(xypos);     %Only keep xypos desired to run
% 
% bkmd.obj.Desc       =   'Apo_20x';  bkmd.obj.Mag        =   20;
% bkmd.obj.WkDist     =   1;          bkmd.obj.RefIndex   =   1;
% bkmd.obj.NA         =   0.75;
% 
% bkmd.cam.PixSizeX   = 	0.65;       bkmd.cam.PixSizeY   =   0.65;  
% bkmd.cam.PixNumX    =   1280;       bkmd.cam.PixNumY    =   1080;
% bkmd.cam.BinSizeX   =   2;          bkmd.cam.BinSizeY   =   2;
% 
% bkmd.exp.Channel    =   {'CFP', 'YFP', 'RFP'};
% bkmd.exp.Exposure   =   [800, 500, 800];   	%Exposure times (ms)
% bkmd.exp.ExVolt     =   [15, 15, 30];       %Percent voltages used
% bkmd.exp.Filter     =   {'Filter_CFP', 'Filter_YFP', 'Filter_Cherry2'};
% bkmd.exp.FPhore     =   {'mTurq2', 'YPet', 'mCherry'};
%  ---  FOR TESTING USE  ---  UNCOMMENT TO TEST PROCEDURE  ---


%% Add any paths necessary for this run
ntl = '\\mcb.ucdavis.edu\Shared\Data\mcb_labs\albeck';
addpath([ntl,'\Code\Image Processing;',ntl,'\bfmatlab']);
%u-Track software
utl = [ntl,'\u-track_2.1.0\software']; 
addpath(utl, [utl,'\mex'], [utl,'\kdtree']);%, [utl,'\bioformats']);   

%Paths needed on workers (Particularly Java ClassPaths)
bfmpath = [ntl,'\bfmatlab\'];
pth.jc = {[bfmpath,'bioformats_package.jar']};
javaclasspath(pth.jc);  %Must include here to have available
%NOTE: To be rid of the "log4j:WARN No appenders foundfor logger" warning,
%   edit classpath.txt, and add the following path to the end of it
%   \\mcb.ucdavis.edu\Shared\Data\mcb_labs\albeck\bfmatlab

       
%% RUN Image Processing
%   Pack parameters to pass to processing function
p = struct('pth',       pth,  'fname',     fname,   'sname',    sname, ...
           'xypos',     xypos,   'nW',       nW, ...
           'op',        op,     'bval',      bval,    'bkmd',     bkmd,   ...
           'bkenforce', {bkenforce},        'bk',     bk, ...
           'xtalk',     xtalk,   'xyshift',  xyshift, 'tsamp',    tsamp, ...
           'disp',      sampdisp);
%   Call processing function
[r, ri, GMD, xtr] = iman_celltrace_parbf2(p, r, ri);
