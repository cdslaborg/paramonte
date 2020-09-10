%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%   MIT License
%%%%
%%%%   ParaMonte: plain powerful parallel Monte Carlo library.
%%%%
%%%%   Copyright (C) 2012-present, The Computational Data Science Lab
%%%%
%%%%   This file is part of the ParaMonte library.
%%%%
%%%%   Permission is hereby granted, free of charge, to any person obtaining a 
%%%%   copy of this software and associated documentation files (the "Software"), 
%%%%   to deal in the Software without restriction, including without limitation 
%%%%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%%%%   and/or sell copies of the Software, and to permit persons to whom the 
%%%%   Software is furnished to do so, subject to the following conditions:
%%%%
%%%%   The above copyright notice and this permission notice shall be 
%%%%   included in all copies or substantial portions of the Software.
%%%%
%%%%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%%%%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%%%%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%%%%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
%%%%   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
%%%%   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
%%%%   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%
%%%%   ACKNOWLEDGMENT
%%%%
%%%%   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
%%%%   As per the ParaMonte library license agreement terms, if you use any parts of 
%%%%   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
%%%%   work (education/research/industry/development/...) by citing the ParaMonte 
%%%%   library as described on this page:
%%%%
%%%%       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2015, Dr. Zdravko Botev
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The University of New South Wales nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
function [bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
    % fast and accurate state-of-the-art
    % bivariate kernel density estimator
    % with diagonal bandwidth matrix.
    % The kernel is assumed to be Gaussian.
    % The two bandwidth parameters are
    % chosen optimally without ever
    % using/assuming a parametric model for the data or any "rules of thumb".
    % Unlike many other procedures, this one
    % is immune to accuracy failures in the estimation of
    % multimodal densities with widely separated modes (see examples).
    % INPUTS: data - an N by 2 array with continuous data
    %            n - size of the n by n grid over which the density is computed
    %                n has to be a power of 2, otherwise n=2^ceil(log2(n));
    %                the default value is 2^8;
    % MIN_XY,MAX_XY- limits of the bounding box over which the density is computed;
    %                the format is:
    %                MIN_XY=[lower_Xlim,lower_Ylim]
    %                MAX_XY=[upper_Xlim,upper_Ylim].
    %                The dafault limits are computed as:
    %                MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
    %                MAX_XY=MAX+Range/4; MIN_XY=MIN-Range/4;
    % OUTPUT: bandwidth - a row vector with the two optimal
    %                     bandwidths for a bivaroate Gaussian kernel;
    %                     the format is:
    %                     bandwidth=[bandwidth_X, bandwidth_Y];
    %          density  - an n by n matrix containing the density values over the n by n grid;
    %                     density is not computed unless the function is asked for such an output;
    %              X,Y  - the meshgrid over which the variable "density" has been computed;
    %                     the intended usage is as follows:
    %                     surf(X,Y,density)
    % Example (simple Gaussian mixture)
    % clear all
    %   % generate a Gaussian mixture with distant modes
    %   data=[randn(500,2);
    %       randn(500,1)+3.5, randn(500,1);];
    %   % call the routine
    %     [bandwidth,density,X,Y]=kde2d(data);
    %   % plot the data and the density estimate
    %     contour3(X,Y,density,50), hold on
    %     plot(data(:,1),data(:,2),'r.','MarkerSize',5)
    %
    % Example (Gaussian mixture with distant modes):
    %
    % clear all
    %  % generate a Gaussian mixture with distant modes
    %  data=[randn(100,1), randn(100,1)/4;
    %      randn(100,1)+18, randn(100,1);
    %      randn(100,1)+15, randn(100,1)/2-18;];
    %  % call the routine
    %    [bandwidth,density,X,Y]=kde2d(data);
    %  % plot the data and the density estimate
    %  surf(X,Y,density,'LineStyle','none'), view([0,60])
    %  colormap hot, hold on, alpha(.8)
    %  set(gca, 'color', 'blue');
    %  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
    %
    % Example (Sinusoidal density):
    %
    % clear all
    %   X=rand(1000,1); Y=sin(X*10*pi)+randn(size(X))/3; data=[X,Y];
    %  % apply routine
    %  [bandwidth,density,X,Y]=kde2d(data);
    %  % plot the data and the density estimate
    %  surf(X,Y,density,'LineStyle','none'), view([0,70])
    %  colormap hot, hold on, alpha(.8)
    %  set(gca, 'color', 'blue');
    %  plot(data(:,1),data(:,2),'w.','MarkerSize',5)
    %
    %  Reference:
    % Kernel density estimation via diffusion
    % Z. I. Botev, J. F. Grotowski, and D. P. Kroese (2010)
    % Annals of Statistics, Volume 38, Number 5, pages 2916-2957.
    global N A2 I
    if nargin<2
        n=2^8;
    end
    n=2^ceil(log2(n)); % round up n to the next power of 2;
    N=size(data,1);
    if nargin<3
        MAX=max(data,[],1); MIN=min(data,[],1); Range=MAX-MIN;
        MAX_XY=MAX+Range/2; MIN_XY=MIN-Range/2;
    end
    scaling=MAX_XY-MIN_XY;
    if N<=size(data,2)
        error('data has to be an N by 2 array where each row represents a two dimensional observation')
    end
    transformed_data=(data-repmat(MIN_XY,N,1))./repmat(scaling,N,1);
    %bin the data uniformly using regular grid;
    initial_data=ndhist(transformed_data,n);
    % discrete cosine transform of initial data
    a= dct2d(initial_data);
    % now compute the optimal bandwidth^2
      I=(0:n-1).^2; A2=a.^2;
     t_star=root(@(t)(t-evolve(t)),N);
    p_02=func([0,2],t_star);p_20=func([2,0],t_star); p_11=func([1,1],t_star);
    t_y=(p_02^(3/4)/(4*pi*N*p_20^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
    t_x=(p_20^(3/4)/(4*pi*N*p_02^(3/4)*(p_11+sqrt(p_20*p_02))))^(1/3);
    % smooth the discrete cosine transform of initial data using t_star
    a_t=exp(-(0:n-1)'.^2*pi^2*t_x/2)*exp(-(0:n-1).^2*pi^2*t_y/2).*a; 
    % now apply the inverse discrete cosine transform
    if nargout>1
        density=idct2d(a_t)*(numel(a_t)/prod(scaling));
        density(density<0)=eps; % remove any negative density values
        [X,Y]=meshgrid(MIN_XY(1):scaling(1)/(n-1):MAX_XY(1),MIN_XY(2):scaling(2)/(n-1):MAX_XY(2));
    end
    bandwidth=sqrt([t_x,t_y]).*scaling; 
end

%#######################################

function [out,time]=evolve(t)
    global N
    Sum_func = func([0,2],t) + func([2,0],t) + 2*func([1,1],t);
    time=(2*pi*N*Sum_func)^(-1/3);
    out=(t-time)/time;
end

%#######################################

function out=func(s,t)
    global N
    if sum(s)<=4
        Sum_func=func([s(1)+1,s(2)],t)+func([s(1),s(2)+1],t); const=(1+1/2^(sum(s)+1))/3;
        time=(-2*const*K(s(1))*K(s(2))/N/Sum_func)^(1/(2+sum(s)));
        out=psi(s,time);
    else
        out=psi(s,t);
    end
end

%#######################################

function out=psi(s,Time)
    global I A2
    % s is a vector
    w=exp(-I*pi^2*Time).*[1,.5*ones(1,length(I)-1)];
    wx=w.*(I.^s(1));
    wy=w.*(I.^s(2));
    out=(-1)^sum(s)*(wy*A2*wx')*pi^(2*sum(s));
end

%#######################################

function out=K(s)
    out=(-1)^s*prod((1:2:2*s-1))/sqrt(2*pi);
end

%#######################################

function data=dct2d(data)
    % computes the 2 dimensional discrete cosine transform of data
    % data is an nd cube
    [nrows,ncols]= size(data);
    if nrows~=ncols
        error('data is not a square array!')
    end
    % Compute weights to multiply DFT coefficients
    w = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
    weight=w(:,ones(1,ncols));
    data=dct1d(dct1d(data)')';
    function transform1d=dct1d(x)

        % Re-order the elements of the columns of x
        x = [ x(1:2:end,:); x(end:-2:2,:) ];

        % Multiply FFT by weights:
        transform1d = real(weight.* fft(x));
    end
end

%#######################################

function data = idct2d(data)
    % computes the 2 dimensional inverse discrete cosine transform
    [nrows,ncols]=size(data);
    % Compute wieghts
    w = exp(i*(0:nrows-1)*pi/(2*nrows)).';
    weights=w(:,ones(1,ncols));
    data=idct1d(idct1d(data)');
    function out=idct1d(x)
        y = real(ifft(weights.*x));
        out = zeros(nrows,ncols);
        out(1:2:nrows,:) = y(1:nrows/2,:);
        out(2:2:nrows,:) = y(nrows:-1:nrows/2+1,:);
    end
end

%#######################################

function binned_data=ndhist(data,M)
    % this function computes the histogram
    % of an n-dimensional data set;
    % 'data' is nrows by n columns
    % M is the number of bins used in each dimension
    % so that 'binned_data' is a hypercube with
    % size length equal to M;
    [nrows,ncols]=size(data);
    bins=zeros(nrows,ncols);
    for i=1:ncols
        [dum,bins(:,i)] = histc(data(:,i),[0:1/M:1],1);
        bins(:,i) = min(bins(:,i),M);
    end
    % Combine the  vectors of 1D bin counts into a grid of nD bin
    % counts.
    binned_data = accumarray(bins(all(bins>0,2),:),1/nrows,M(ones(1,ncols)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t=root(f,N)
    % try to find smallest root whenever there is more than one
    N=50*(N<=50)+1050*(N>=1050)+N*((N<1050)&(N>50));
    tol=10^-12+0.01*(N-50)/1000;
    flag=0;
    while flag==0
        try
            t=fzero(f,[0,tol]);
            flag=1;
        catch
            tol=min(tol*2,.1); % double search interval
        end
        if tol==.1 % if all else fails
            t=fminbnd(@(x)abs(f(x)),0,.1); flag=1;
        end
    end
end
