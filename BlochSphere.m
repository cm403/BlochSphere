% Bloch Sphere
% for i = 1
while 1
    pause(.1)
    theta = mod(theta +pi/111,2*pi);
    phi = mod(phi+pi/55,2*pi);

    Psi = [cos(theta/2)  exp(1i*phi)*sin(theta/2)]
    norm(Psi)
    thetaphi = ([theta phi]);

    azel = thetaphi(2:-1:1).*[1 -1] + deg2rad([0  90]);

    [psix psiy psiz] = sph2cart(azel(1), azel(2), 1);
    psixyz = [psix psiy psiz];
    %% Draw Sphere
    ph = deg2rad(-90:2:90);
    th = deg2rad(-180:2:180);

    C = ones(181,91);
    %C(:,:,:) = 1;

    th = th;
    ph = ph;
    [ph, th] = meshgrid(ph, th);

    color = [0 0 1];


    Beta = deg2rad(90);
    clf


    for i = 1


        X = C(:,:,i).*sin(th).*cos(ph);
        Y = C(:,:,i).*sin(th).*sin(ph);
        Z = C(:,:,i).*cos(th);

        Xr = X;
        Yr = cos(Beta)*Y + -sin(Beta)*Z;
        Zr = sin(Beta)*Y + cos(Beta)*Z;

        C(:,:,i) = sqrt(Xr.^2 + Yr.^2 +Zr.^2);


        figure(2);
        X = C(:,:,i).*sin(th).*cos(ph);
        Y = C(:,:,i).*sin(th).*sin(ph);
        Z = C(:,:,i).*cos(th);

        mesh(Xr,Yr,Zr,'FaceColor',color(i,:),'EdgeColor','none');
        alpha(.2)
        hold on
    end
    view(120,30)

    p1 = [0 0 0];                         % First Point
    p2 = [1 0 0;
          0 1 0
          0 0 1];                         % Second Point
    basis = {'\midv\rangle';'\midw\rangle';'\mid0\rangle'};
    for i = 1:3
        dp = p2(i,:)-p1;                         % Difference

        quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0)

    %     text(p1(1),p1(2), sprintf('(%.0f,%.0f)',p1))
        text(p2(i,1),p2(i,2),p2(i,3), sprintf('%s',basis{i}))

    %     text(p2(i,1),p2(i,2),p2(i,3), sprintf('(%.0f,%.0f,%.0f)',p2(i,:)))
    end


    xlabel('x')
    ylabel('y')
    zlabel('z')

    %% Plot Psi
    p2 = [psixyz;     
          ] ;                   % Second Point

    basis = {'\Psi',''};
    for i = 1:size(p2,1)
        plot3([p1(1);p2(i,1)],[p1(2);p2(i,2)],[p1(3);p2(i,3)])
        text(p2(i,1),p2(i,2),p2(i,3), sprintf('%s',basis{i}))

    end

    % Plot Psi Components
    p2 = [psixyz.*[1 1 0]
          psixyz
          ] ;                   % Second Point
    p1 = [0 0 0;
          psixyz.*[1 1 0]]; 

    for i = 1:size(p2,1)
        plot3([p1(i,1);p2(i,1)],[p1(i,2);p2(i,2)],[p1(i,3);p2(i,3)],'k--')

    end

    % Plot Psi Angles
    p2 = [psixyz.*[1 1 0]/3
          psixyz/3
          ] ;                   % Second Point
    p1 = [[1 0 0]/3;
          [0 0 1]/3]; 

    s = linspace(0,phi,100);
    t = linspace(0,theta,100);
    if norm(psixyz(1:2)) < 1/3
        r1 = norm(psixyz(1:2));
    else
        r1 = 1/3;
    end
    r2 = 1/3;
    X = [r1*cos(s)' r2*sin(t)'.*cos(-phi)'];
    Y = [r1*sin(s)' r2*sin(t)'.*-sin(-phi)'];
    Z = [0*s' r2*cos(t)'];
    ang = {'\phi','\theta'}
    for i = 1:2
        plot3(X(:,i),Y(:,i),Z(:,i),'k-')
        text(X(50,i)',Y(50,i)',Z(50,i)', sprintf('%s',ang{i}))
    end

    %% Title

    tl1 = sprintf('%s_u = ( %0.2f + %0.2fi ) %s0%s + ( %0.2f + %0.2fi ) %s1%s',...
        '\psi',Psi(1),Psi(1)*1i,'\mid','\rangle',Psi(2),Psi(2)*1i,'\mid','\rangle');
    Psiv = 1/sqrt(2)*[Psi(1)+Psi(2) Psi(1)-Psi(2)];
    tl2 = sprintf('%s_v = ( %0.2f + %0.2fi ) %s0%s + ( %0.2f + %0.2fi ) %s1%s',...
        '\psi',Psiv(1),Psiv(1)*1i,'\mid','\rangle',Psiv(2),Psiv(2)*1i,'\mid','\rangle');
    Psiw = 1/sqrt(2)*[Psi(1)-1i*Psi(2) Psi(1)+1i*Psi(2)];
    tl3 = sprintf('%s_w = ( %0.2f + %0.2fi ) %s0%s + ( %0.2f + %0.2fi ) %s1%s',...
        '\psi',Psiw(1),Psiw(1)*1i,'\mid','\rangle',Psiw(2),Psiw(2)*1i,'\mid','\rangle');

    tl4 = sprintf('%s = %0.1f%s, %s = %0.1f%s','\theta',rad2deg(theta),char(176),'\phi',rad2deg(phi),char(176))

    title(sprintf('%s\n%s\n%s\n%s',tl1,tl2,tl3,tl4))
end



