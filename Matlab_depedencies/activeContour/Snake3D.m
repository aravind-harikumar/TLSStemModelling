function FV=Snake3D(I,FV,Options)
% This function SNAKE implements the basic snake segmentation. A snake is an
% active (moving) contour, in which the points are attracted by edges and
% other boundaries. To keep the contour smooth, an membrame and thin plate
% energy is used as regularization.
%
% OV=Snake3D(I,FV,Options)
%
% inputs,
%   I : An Image of type double preferable ranged [0..1]
%   FV : Structure with triangulated mesh, with list of faces FV.faces N x 3
%        and list of vertices M x 3
%   Options : A struct with all snake options
%
% outputs,
%   OV : Structure with triangulated mesh of the final surface
%
% options (general),
%  Option.Verbose : If true show important images, default false
%  Options.Gamma : Time step, default 1
%  Options.Iterations : Number of iterations, default 100
%
% options (Image Edge Energy / Image force))
%  Options.Sigma1 : Sigma used to calculate image derivatives, default 2
%  Options.Wline : Attraction to lines, if negative to black lines otherwise white
%                    lines , default 0.04
%  Options.Wedge : Attraction to edges, default 2.0
%  Options.Sigma2 : Sigma used to calculate the gradient of the edge energy
%                    image (which gives the image force), default 2
%
% options (Gradient Vector Flow)
%  Options.Mu : Trade of between real edge vectors, and noise vectors,
%                default 0.2. (Warning setting this to high >0.5 gives
%                an instable Vector Flow)
%  Options.GIterations : Number of GVF iterations, default 0
%  Options.Sigma3 : Sigma used to calculate the laplacian in GVF, default 1.0
%
% options (Snake)
%  Options.Alpha : Membrame energy  (first order), default 0.2
%  Options.Beta : Thin plate energy (second order), default 0.2
%  Options.Delta : Baloon force, default 0.1
%  Options.Kappa : Weight of external image force, default 2
%  Options.Lambda : Weight which changes the direction of the image
%                   potential force to the direction of the surface
%                   normal, value default 0.8 (range [0..1])
%                   (Keeps the surface from self intersecting)
%
% Literature:
%   - Michael Kass, Andrew Witkin and Demetri TerzoPoulos "Snakes : Active
%       Contour Models", 1987
%   - Christoph Lurig, Leif Kobbelt, Thomas Ertl, "Hierachical solutions
%       for the Deformable Surface Problem in Visualization"


% Process inputs
defaultoptions=struct('Verbose',false,'Wline',0.04,'Wedge',2,'Sigma1',2,'Sigma2',2,'Alpha',0.2,'Beta',0.2,'Delta',0.1,'Gamma',1,'Kappa',2,'Iterations',100,'GIterations',0,'Mu',0.2,'Sigma3',1,'Lambda',0.8);
if(~exist('Options','var')),
    Options=defaultoptions;
else
    tags = fieldnames(defaultoptions);
    for i=1:length(tags)
        if(~isfield(Options,tags{i})), Options.(tags{i})=defaultoptions.(tags{i}); end
    end
    if(length(tags)~=length(fieldnames(Options))),
        warning('snake:unknownoption','unknown options found');
    end
end

% Convert input to single if xintxx
if(~strcmpi(class(I),'single')&&~strcmpi(class(I),'double'))
    I=single(I);
end

% The surface faces must always be clockwise (because of the balloon force)
FV=MakeContourClockwise3D(FV);

% Transform the Image into an External Energy Image
Eext = ExternalForceImage3D(I,Options.Wline, Options.Wedge,Options.Sigma1);

% Make the external force (flow) field.
Fx=ImageDerivatives3D(Eext,Options.Sigma2,'x');
Fy=ImageDerivatives3D(Eext,Options.Sigma2,'y');
Fz=ImageDerivatives3D(Eext,Options.Sigma2,'z');

Fext(:,:,:,1)=-Fx*2*Options.Sigma2^2;
Fext(:,:,:,2)=-Fy*2*Options.Sigma2^2;
Fext(:,:,:,3)=-Fz*2*Options.Sigma2^2;

% Do Gradient vector flow, optimalization
Fext=GVFOptimizeImageForces3D(Fext, Options.Mu, Options.GIterations, Options.Sigma3);

% Show the image, contour and force field
% if(Options.Verbose)
%     drawnow; pause(0.1);
    %      h=figure; set(h,'render','opengl')
    %      subplot(2,3,1),imshow(squeeze(Eext(:,:,round(end/2))),[]);
    %      subplot(2,3,2),imshow(squeeze(Eext(:,round(end/2),:)),[]);
    %      subplot(2,3,3),imshow(squeeze(Eext(round(end/2),:,:)),[]);
    %      subplot(2,3,4),imshow(squeeze(Fext(:,:,round(end/2),:))+0.5);
    %      subplot(2,3,5),imshow(squeeze(Fext(:,round(end/2),:,:))+0.5);
    %      subplot(2,3,6),imshow(squeeze(Fext(round(end/2),:,:,:))+0.5);
    %      h=figure; set(h,'render','opengl'); hold on;
    %patch(i,'facecolor',[1 0 0],'facealpha',0.1);
    ind=find(I(:)>0);
    [ix,iy,iz]=ind2sub(size(Eext),ind);
%     subplot(2,4,2);
%     plot3(ix,iy,iz,'r.','MarkerSize',3);
%     subplot(1,3,3);
%     plot3(ix,iy,iz,'r.','MarkerSize',3);
    %      camproj perspective; view([-45 15]);
    %      set(gcf, 'Position', get(0, 'Screensize'));
%     pause(0.5);
%     hold on;
%     h=patch(FV,'facecolor',[0 0 0.5],'facealpha',0.1);
%     drawnow;
%     pause(0.2);
% end
% axis([0 150 0 150 0 400]);
% camproj perspective; %view([-45 25]);
% pbaspect([2 2 10])
% view([-90 0]);
%      axis equal
%      axis([0 25 0 25 225 555]);
% Make the interal force matrix, which constrains the moving points to a
% smooth contour
%  figure;
S=SnakeInternalForceMatrix3D(FV,Options.Alpha,Options.Beta,Options.Gamma);
set(gcf, 'Position', get(0, 'Screensize'));
for i=1:Options.Iterations
    FV=SnakeMoveIteration3D(S,FV,Fext,Options.Gamma,Options.Kappa,Options.Delta,Options.Lambda);
    % Show current contour
%     if(Options.Verbose)
%         if(ishandle(h))
%             delete(h);
%             FV.vertices = real(FV.vertices);
%             subplot(2,4,2);
%             plot3(ix,iy,iz,'r.','MarkerSize',3);
%             hold on;
%             h=patch(FV,'facecolor',[1 0 0],'facealpha',0.1);
%             drawnow; 
%             axis equal;
%             axis([0 200 0 200 10 20]);
%             camproj perspective; view([90 90]);
%             hold off;
%             figure;
            subplot(2,4,2);
            plot3(ix,iy,iz,'r.','MarkerSize',3);
            hold on;
            h=patch(FV,'facecolor',[1 0 0],'facealpha',0.1);
            drawnow; 
            axis equal;
            axis([0 size(I,1) 0 size(I,2) 85 145]);
            camproj perspective; view([90 90]);
            title(strcat( 'Cross Section: ', num2str('1.25 to 1.35 m')))
            pause(0.1);
            hold off;
             
            subplot(2,4,3);
            plot3(ix,iy,iz,'r.','MarkerSize',3);
            hold on;
            h=patch(FV,'facecolor',[1 0 0],'facealpha',0.1);
            drawnow; 
            axis equal;
            axis([0 size(I,1) 0 size(I,2) 145 185]);
            camproj perspective; view([90 90]);
            title(strcat( 'Cross Section: ', num2str('1.45 to 1.55 m')))
            pause(0.1);
            hold off;

            subplot(2,4,4);
            plot3(ix,iy,iz,'r.','MarkerSize',3);
            hold on;
            h=patch(FV,'facecolor',[1 0 0],'facealpha',0.1);
            drawnow; 
            axis equal;
            axis([0 size(I,1) 0 size(I,2) 1 250]);
            camproj perspective; view([45 45]);
            title(strcat( 'Cross Section: ', num2str('0.00 to 2.55 m')))
            pause(0.1);
            hold off;             
            %              axis([0 25 0 25 225 555]);
            %             pause(0.5);
%         end
%     end
end


