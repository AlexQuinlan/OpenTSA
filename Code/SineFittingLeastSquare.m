function [a,p,o,e] = SineFittingLeastSquare(f,t,d,pf,ff)
            % Input:
            %   f: reference frequency
            %   t: time base
            %   d: input matrix 3D
            % Output:
            %   a: amplitude
            %   p: phase
            %   o: offset
            %   e: RMSD error
            t_ = reshape(t,[],1);
            A = [cos(2*pi*f*t_) -sin(2*pi*f*t_) ones(size(t_))];
            a = zeros(size(d,1),size(d,2));
            p = zeros(size(d,1),size(d,2));
            o = zeros(size(d,1),size(d,2));
            e = zeros(size(d,1),size(d,2));
            
            if length(size(d)) > 2 
                for x = 1:size(d,2)
                    for y = 1:size(d,1)
                        d1 = permute(d(y,x,:),[3 1 2]);
                        d2 = pf(d1);
                        d3 = ff(d2);        
%                         d3 = d1;
                        te = A\d3;                    
                        a(y,x) = sqrt(te(1)^2+te(2)^2);
                        p(y,x) = atan2(te(2),te(1));
                        o(y,x) = te(3);
                        e(y,x) = sqrt( mean( (d3-A*te).^2 ) );
    %                         [a(y,x) p(y,x) o(y,x) e(y,x)]
                    end
                end
            else
                                    
                te = A\d;                    
                a = sqrt(te(1)^2+te(2)^2);
                p = atan2(te(2),te(1));
                o = te(3);
                e = sqrt( mean( (d-A*te).^2 ) ); 
            end
        end