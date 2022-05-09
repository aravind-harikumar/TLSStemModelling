function [y,x] = plotellipise(pmval,ploton)
    
%     % S = solve(x^(5/2) == 8^(sym(10/3))) returns all three complex solutions:
%     arraRoot = [];
%     for y =1:1:100
%     %syms x
%     %opp = solve(pmval(1)*x^2 + pmval(2)*x*y + pmval(3)*y^2 + pmval(4)*x + pmval(5)*y + pmval(6) == 0,'IgnoreAnalyticConstraints', true);
%     cc = pmval(3)*y^2 +  pmval(5)*y + pmval(6);    
%     bb = (pmval(2)*y + pmval(4));    
%     aa = pmval(1);    
%     dd = sqrt(bb^2 - 4*aa*cc);
%     rootss(1) = real(( -bb + dd ) / (2*aa)); imrootss(1) = imag(( -bb + dd ) / (2*aa));
%     rootss(2) = real (( -bb - dd ) / (2*aa)); imrootss(2) = imag (( -bb - dd ) / (2*aa));
%     
%     if(rootss(1) ~= rootss(2))    
%     arraRoot = [arraRoot; rootss(1) y ; rootss(2) y];
%     end
%     end
%     
%     plot(arraRoot(:,2), arraRoot(:,1),'*');
%     
%     
%    hold on;
    
   A = pmval(1); B = pmval(2); C = pmval(3); D = pmval(4); E = pmval(5); F= pmval(6);
   e = 4*A*C-B^2; if e<=0, error('This conic is not an ellipse.'), end
   x0 = (B*E-2*C*D)/e; y0 = (B*D-2*A*E)/e;   % Ellipse center
   F0 = -2*(A*x0^2+B*x0*y0+C*y0^2+D*x0+E*y0+F);
   g = sqrt((A-C)^2+B^2); a = F0/(A+C+g); b = F0/(A+C-g);
   if (a<=0)|(b<=0), error('This is a degenerate ellipse.'), end
   a = sqrt(a);  b = sqrt(b); % Major & minor axes
   t = 1/2*atan2(B,A-C); ct = cos(t); st = sin(t);   % Rotation angle
   p = linspace(0,2*pi,100); cp = cos(p); sp = sin(p);   % Variable parameter
   x = x0+a*ct*cp-b*st*sp; y = y0+a*st*cp+b*ct*sp;   % Generate points on ellipse
   if(ploton)
   plot(y,x,'r-', 'MarkerSize', 15), axis equal   % Plot them
   end
   x = x'; y = y';
    
%     y = 1:1:200;
%     temp = -(pmval(2)*x*y + pmval(3)*y^2 + pmval(4)*x + pmval(5)*y + pmval(6));
%     op = sqrt(temp/pmval(1));
end