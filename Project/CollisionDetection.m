function [n_obj,mass,x,y,vx,vy,px,py,rad] = CollisionDetection(n_obj,mass,x,y,vx,vy,px,py,rad)
    
pairr=[];
pir=[];
k=0;
    for ii = 1:numel(mass)-1
        for jj=ii+1:numel(mass)
            if(length(find(ii==pairr(:)))>0)
                continue;
            end
            if(length(find(jj==pairr(:)))>0)
                continue;
            end

            dist=sqrt((x(ii)-x(jj)).^2+(y(ii)-y(jj)).^2);
            radii=rad(ii)+rad(jj);
            idx=find(dist<radii);

            if(length(idx)>0)
                k=k+1;
                pir(k,:)=[ii jj];
                pairr=[pairr ii jj];
                cflg=1;
                break;
            end
        end
    end

    for ik = 1:size(pir,1)
        ii=pir(ik,1);
        jj=pir(ik,2);
        mass(end+1)=mass(ii)+mass(jj);
        edge_a=rad(ii)/rad(jj);
        edge_b=rad(jj)/rad(ii);
        
        if(edge_a>=2)
            x(end+1)=x(ii);
            y(end+1)=y(ii);
        elseif(edge_b>=2)
            x(end+1)=x(jj);
            y(end+1)=y(jj);
        else
            x(end+1)=(x(ii)+x(jj))/2;
            y(end+1)=(y(ii)+y(jj))/2;
        end
        vx(end+1)=(px(ii)+px(jj))/mass(end);
        vy(end+1)=(py(ii)+py(jj))/mass(end); 
        px(end+1)=(mass(end)*vx(end));
        py(end+1)=(mass(end)*vy(end));
    end

    mass(pairr)=[];
    x(pairr)=[];
    y(pairr)=[];
    vx(pairr)=[];
    vy(pairr)=[];
    px(pairr)=[];
    py(pairr)=[];

    rad= sqrt(mass);
end