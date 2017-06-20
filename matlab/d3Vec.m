classdef d3Vec

	properties
		x
		y
		z
	end

	methods

		function obj = d3Vec(varargin)
			%-- The constructor for the d3Vec class. With no argument it makes the
			%fields empty. With one input, if the class of the input is d3Vec, it makes
			%a copy of input, if the input is a vector/matrix it is assumed that the
			%coordinate points are ordered in [X;Y;Z] or [X Y Z] order.
			%Three inputs are used to assign each input to a field of the object.
			switch nargin
				case 0
					X = []; Y = []; Z = [];
				case 1
					a = varargin{1};
					if(isa(a,'d3Vec'))
						X = a.x; Y = a.y; Z = a.z;
					else
						a = reshape(a,[],3);
						X = a(:,1); Y = a(:,2); Z = a(:,3);
					end
				case 3
					X = varargin{1};
					Y = varargin{2};
					Z = varargin{3};
			end
			obj.x = X;
			obj.y = Y;
			obj.z = Z;

		end

		function vec = vecForm(obj)
                  for ii=1:length(obj)
                    vec(:,ii) = [obj(ii).x;obj(ii).y;obj(ii).z];
                  end
		end
		
		function val = mtimes(mat,vec)
			%-- scalar multiplication.
			vec.x = mat*vec.x;
			vec.y = mat*vec.y;
			vec.z = mat*vec.z;
			val = vec;
		end

		function val = mrdivide(vec,mat)
			%-- scalar division
			vec.x = vec.x/mat;
			vec.y = vec.y/mat;
			vec.z = vec.z/mat;
			val = vec;
		end

		function vec1 = plus(vec1,vec2)
			vec1.x = vec1.x+vec2.x;
			vec1.y = vec1.y+vec2.y;
			vec1.z = vec1.z+vec2.z;
		end

		function vec1 = minus(vec1,vec2)
			vec1.x = vec1.x-vec2.x;
			vec1.y = vec1.y-vec2.y;
			vec1.z = vec1.z-vec2.z;
		end

		function val = cross(vec1,vec2)
			val = d3Vec;
			val.x =  vec1.y.*vec2.z - vec1.z.*vec2.y;
			val.y = -vec1.x.*vec2.z + vec1.z.*vec2.x;
			val.z =  vec1.x.*vec2.y - vec1.y.*vec2.x;
		end

		function val = dot(vec1,vec2)
			val = vec1.x.*conj(vec2.x) + vec1.y.*conj(vec2.y) + vec1.z.*conj(vec2.z);
		end

		function vec = times(scalar,vec)
			%-- elementwise multioication .*
			if(isa(scalar,'d3Vec'))
				vec.x = scalar.x.*vec.x;
				vec.y = scalar.y.*vec.y;
				vec.z = scalar.z.*vec.z;
			else
				vec.x = scalar.*vec.x;
				vec.y = scalar.*vec.y;
				vec.z = scalar.*vec.z;
			end
		end

		function vec = rdivide(vec,scalar)
			%-- elementwise division ./
			if(isa(scalar,'d3Vec'))
				vec.x = vec.x./scalar.x;
				vec.y = vec.y./scalar.y;
				vec.z = vec.z./scalar.z;
			else
				vec.x = vec.x./scalar;
				vec.y = vec.y./scalar;
				vec.z = vec.z./scalar;
			end
		end

		function vec = ctranspose(vec)
			vec.x = vec.x';
			vec.y = vec.y';
			vec.z = vec.z';
		end

		function vec = transpose(vec)
			vec.x = vec.x.';
			vec.y = vec.y.';
			vec.z = vec.z.';
		end

		function res = sum(vec,dim)
			if(nargin<2), dim = 1; end
				
			res = d3Vec;
			if(dim==1)
				res.x = sum(vec.x);
				res.y = sum(vec.y);
				res.z = sum(vec.z);
			else
				res = vec.x+vec.y+vec.z;
			end
			
		end
		
		function vec = abs(vec)
			vec.x = abs(vec.x);
			vec.y = abs(vec.y);
			vec.z = abs(vec.z);
		end
		
		function vec = shAna(vec)
			for ii=1:length(vec)
                          c = shAna([vec(ii).x vec(ii).y vec(ii).z]);
                          vec(ii).x = c(:,1);
                          vec(ii).y = c(:,2);
                          vec(ii).z = c(:,3);
                        end
		end
		
		function vec = shSyn(vec,isReal)
			if(nargin<2), isReal = [];end
			for ii=1:length(vec)
                          c = shSyn([vec(ii).x vec(ii).y vec(ii).z],isReal);
                          vec(ii).x = c(:,1);
                          vec(ii).y = c(:,2);
                          vec(ii).z = c(:,3);
                        end
		end
                		
		function vec = interpsh(vec,newSize)
			c = interpsh([vec.x vec.y vec.z],newSize);
			vec.x = c(:,1);
			vec.y = c(:,2);
			vec.z = c(:,3);
		end

		function [vec mat] = filterSh(vec,filtFreq)
			if(nargout>1)
                          [c mat] = filterSh([vec.x vec.y vec.z],filtFreq);
                        else
                          c = filterSh([vec.x vec.y vec.z],filtFreq);
                        end  
			vec.x = c(:,1);
			vec.y = c(:,2);
			vec.z = c(:,3);
		end
		
		function vec = DmSH(vec,dv,du)
			c = DmSH([vec.x vec.y vec.z],dv,du);
			vec.x = c(:,1);
			vec.y = c(:,2);
			vec.z = c(:,3);
		end

                function obj=reshape(obj,nx,ny)
                        obj.x = reshape(obj.x,nx,ny);
                        obj.y = reshape(obj.y,nx,ny);
                        obj.z = reshape(obj.z,nx,ny);
                end
		
		function quiver3(v1,v2,varargin)
			quiver3(v1.x,v1.y,v1.z,v2.x,v2.y,v2.z,varargin{:});
		end
		
		function nv = norm(vec)
			nv = sqrt(dot(vec,vec));
		end                
                
                function obj = meshgrid(obj, varargin)
                  [obj.x obj.y obj.z] = meshgrid(varargin{:});
                end

                function obj = movePole(obj, varargin)
                  cart = movePole(reshape(obj.vecForm,[],3), varargin{:});
                  obj.x = cart(:,1);
                  obj.y = cart(:,2);
                  obj.z = cart(:,3);
                  
                end
	end
end
