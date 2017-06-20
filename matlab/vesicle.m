classdef vesicle < handle 
    
    properties
        kappa;    %the bending modulus
        tension;  %the array holding tension
    end
    
    properties(SetObservable) %i.e. changes of these will cause other
                              %memebers to get updated

        cart = d3Vec;        %Cartesian coordinate of the points. cart is of
                             %type d3Vec that is user defined.
        
        upFreqCoeff;     %upsampling rate for differentiation
        filterFreqCoeff; %filtering rate for differentiation
        viscCont;        %viscosity contrast
        stokesOpStale = true;   %flag indicating whether the precomputed
                                %stokes is stale or not (in case the stokes
                                %operator is stored as a matrix)
        dblLayerOpStale = true; %same as above for double layer
    end
    
    properties(SetObservable, SetAccess = 'protected')
        p             %The number of spherical harmonics frequencies. The
                      %number of points is 2*p*(p+1).
    end
    
    properties(SetAccess = 'protected')
        shc = d3Vec;         %Spherical harmonic coefficients of Cartesian
                             %coordinates. shc is of type d3Vec (user defined
                             %type).
        
        geoProp;             %Geometric properties of the surface. It
                             %contains first and second derivative, first and
                             %second fundamental coefficient, normal vector,
                             %Gaussian and mean curvature.
        
        upFreq;              % upsampling frequency for differentiation  = upFreqCoeff * p
        filterFreq;          % filtering frequency for differentiation = filterFreqCoeff * p
        velCoeff;            % 2/(1+viscCont)
        dblLayerCoeff = 0;   % (1-viscCont)
    end
        
    methods

        %Constructor for the objects of type vesicle. It can be initialized
        %by giving the Cartesian coordinate of the points. Coordinate points
        %are stored in (and can be initialized with) an object of type d3Vec.
        function obj=vesicle(cartInit, kappaIn, upFreqCoeffIn, ...
                filterFreqCoeffIn, viscContIn)
            
            if(isa(cartInit,'vesicle'))
                obj.kappa = cartInit.kappa;
                obj.upFreqCoeff = cartInit.upFreqCoeff;
                obj.filterFreqCoeff = cartInit.filterFreqCoeff;
                obj.viscCont = cartInit.viscCont;
                obj.cart = cartInit.cart;
            else
                if(nargin>1)
                    obj.kappa = kappaIn;
                else
                    obj.kappa = 1;
                end
                if(nargin>2)
                    obj.upFreqCoeff = upFreqCoeffIn;
                else
                    obj.upFreqCoeff = 1;
                end
                if(nargin>3)
                    obj.filterFreqCoeff = filterFreqCoeffIn;
                else
                    obj.filterFreqCoeff = 1;
                end
            
                if(nargin>4)
                    obj.viscCont = viscContIn;
                else
                    obj.viscCont = 1;
                end
                obj.cart = d3Vec(cartInit);
            end
        end
        
        function obj=set.cart(obj,val)
            %Set method for cart field. cart field is of type d3Vec and can
            %be set to any object of this type.
            
            if(isa(val,'d3Vec'))
                obj.cart = val;
            else
                error(['Surface''s Cartesian coordinates can be set by an object ' ...
                    'of class d3Vec.'])
            end
       
            % Updating the rest of the params
            obj.p = (sqrt(2*size(obj.cart.x,1)+1)-1)/2;
            obj.shc = obj.cart.shAna();
            obj.geoProp = calcGeoProp(obj);
            obj.stokesOpStale = true;
            obj.dblLayerOpStale = true;
            if(~isempty(obj.tension))
                obj.tension = interpsh(obj.tension, obj.p);
            else
                obj.tension = zeros(size(obj.cart.x));
            end
        end
        
        function obj=set.upFreqCoeff(obj,val)
            if( size( val ) ~= size( obj ) )
              val = repmat( val, size( obj ) );
            end
            for ii=1:length(obj)
              flag = (obj(ii).upFreqCoeff ~= val(ii));
              obj(ii).upFreqCoeff = val(ii); 
              if(flag), obj(ii).resample; end
            end
        end
        
        function obj=set.filterFreqCoeff(obj,val)
            if( size( val ) ~= size( obj ) )
              val = repmat( val, size( obj ) );
            end
            for ii=1:length(obj)
              flag = (obj(ii).filterFreqCoeff ~=val(ii));
              obj(ii).filterFreqCoeff = val(ii);
              if(flag), obj(ii).resample; end
            end
        end
        
        function obj=set.viscCont(obj, viscCont)
            for ii=1:length(obj)
              obj(ii).viscCont = viscCont;
              obj(ii).dblLayerCoeff = (1-viscCont);
              obj(ii).velCoeff = 2/(1+viscCont);
            end
        end
        
        function f = get.upFreq(obj)
            f = ceil(obj.p*obj.upFreqCoeff);
        end
        
        function f = get.filterFreq(obj)
            f = floor(obj.p*obj.filterFreqCoeff);
        end
        
        function [Sf callsRet] = stokesMatVec(obj, f, varargin)
       % Stokes matvec wrapper     
            persistent calls
            
            if ( isempty(calls) )
                calls.n = 0;
                calls.t = 0;
            end
            
            if(nargin==1)
                Sf = [];
                callsRet = calls;
                return;
            end
            
            tt = clock;
            Sf = kernelS(f, obj, varargin{:});

            calls.t = calls.t + etime(clock,tt);
            if(~isempty(f)), calls.n = calls.n + 1;end
            if(nargout>1) callsRet = calls; end
        end

        function Sshc = stokesMatVecSh(obj, shc, varargin)
            %shc to shc map for stokes matvec
            Sf = obj.stokesMatVec(shSyn(shc, false), varargin{:});
            Sshc = shAna(Sf);
        end
        
        function [Df numcalls] = doubleLayerMatVec(obj, f, varargin)
            %Double layer operator wrapper
              
            persistent ncalls;
            if ( isempty(ncalls) ), ncalls = 0; end
            
            ncalls = ncalls + 1;
            numcalls = ncalls;
            
            if ( ~any( [ obj.viscCont ] ~= 1 ) )
                for ii=1:length(obj)
                  Df(ii) = 0 * f(ii);
                end
                return;
            end
            
            Df = kernelD(f, obj, varargin{:});
            for ii=1:length(obj)
              Df(ii) = obj(ii).dblLayerCoeff * Df(ii);
            end
        end
                
        function Dshc = doubleLayerMatVecSh(obj, shc, varargin)
            % shc to shc map for double layer
              Df = obj.doubleLayerMatVec(shSyn(shc, false), varargin{:});
            Dshc = shAna(Df);
        end
        
        function vel = tensionStokesSh(obj, sig, isReal, varargin)
            % The tension operator L (shc to shc)
            if(nargin<3), isReal = false; end;
            vel = obj.stokesMatVec(obj.geoProp.tensionOp(shSyn(sig, isReal)), ...
                                   varargin{:});
        end
        
        function div = divTensionStokesSh(obj, sig, isReal, varargin)
            if(nargin<3), isReal = false; end 
            div = shAna(obj.geoProp.Div(obj.tensionStokesSh(sig, isReal, varargin{:})));
        end
        
        function hOut = plot(obj,varargin)
            h=plotb(obj.getAllCart(),varargin{:}); drawnow;
            if(nargout>0), hOut = h; end
        end
        
        function obj = resample(obj,newFreq)
            if(nargin<2), newFreq = obj.p; end
            if(~isempty(obj.cart))
                obj.cart = obj.cart.interpsh(newFreq);
            end
        end
        
        function X = movePole(obj, thetaIdx, phiIdx)
            XX = movePoleShcMatrix(reshape(obj.cart.vecForm, [], 3),...
                thetaIdx, phiIdx);
            for ii=1:length(XX)
                X(ii) = d3Vec(XX{ii});
            end
        end
        
        function obj = filter(obj,filtFreq)
            obj.cart = obj.cart.filterSh(filtFreq);
        end
        
        function v = area(obj)
            [a v] = reducedVolume(obj);
        end
        
        function a = volume(obj)
            a = reducedVolume(obj);
        end
        
        function cent = centerOfMass(obj)
            cent = getCenter(obj);
        end
        
        function ra = reducedVol(obj)
            [v a ra] = reducedVolume(obj);
        end
        
        function II = momentOfInertiaVolume(obj)
            C = obj.centerOfMass();
            R = obj.cart-C;
            S = .5*dot(R,R);
            w = .5*R.*R.*R;
            nor = obj.geoProp.nor;
            II = eye(3);
            for i=1:3
                ei=d3Vec(II(:,i));                               
                f = integrateOverS(obj, dot(w,nor).*ei - dot(R,ei).*S.*nor);
                II(:,i) = f.vecForm;
            end            
            II = II/obj.volume(); %normalize with volume (equivalent to whole mass)
            %expect .4 for sphere and [1 1 .4] for ellipse (1,1,2)
        end
        
        function II = momentOfInertiaShell(obj)
            C = obj.centerOfMass();
            R = obj.cart-C;
            a = dot(R,R);
            II = eye(3);
            for i=1:3
                ei=d3Vec(II(:,i)); 
                f = integrateOverS(obj, a*ei-dot(R,ei).*R);
                II(:,i) = f.vecForm;
            end            
            II = II/obj.area(); %normalize with area (equivalent to shell mass)
        end
        
        function t = tau(obj)
            t = (obj.area/4/pi)^1.5/obj.kappa;
        end
        
        function X = getAllCart(obj)
            for ii=1:length(obj)
                X(ii) = obj(ii).cart;
            end
        end
        
        function ten = getAllTension(obj)
            for ii=1:length(obj)
                ten{1,ii} = obj(ii).tension;
            end
        end
        
        function xs = getAllCartAndTensionSh(obj)
            nv = length(obj);
          for ii=1:nv
                xs(:,ii) = [obj(ii).cart.vecForm; obj(ii).tension];
            end
            xs = reshape(xs, [], 4*nv);
            xs = shAna(xs);
            xs = reshape(xs, [], nv);
        end
        
        function obj = setAllCartAndTensionSh(obj, xs)
            nv = length(obj);
            xs = reshape(xs, [], 4*nv);
            xs = shSyn(xs, true);
            for ii=0:nv-1
                obj(ii+1).cart = d3Vec(xs(:,4*ii+1:4*ii+3));
                obj(ii+1).tension = xs(:,4*ii+4);
            end
        end
    end
end
