classdef FixedLengthQueue < handle
    properties
        l
    end
    
    properties (Access=private)
        v
        n
    end
    
    methods
        function this = FixedLengthQueue(l)
            this.v = cell(0);
            this.l = l;
            this.n = 0;
        end
        function add(this, p)
            if this.n == this.l
                this.v = circshift(this.v,-1);
                this.v{end} = p;
            else
                this.v{end+1} = p;
                this.n = this.n + 1;
            end
        end
        
        function p = remove(this)
            if this.n ~= 0
                p = this.v{1};
                if this.n > 1
                    this.v = this.v(2:end);
                else
                    this.v = cell(0);
                end
                this.n = this.n - 1;
            else
                error('Empty queue.')
            end
        end
        
        function p = peek(this)
            if this.n ~= 0
                p = this.v{1};
            else
                error('Empty queue.')
            end
        end
        
        function pVec = peekAll(this)
            if this.n ~= 0
                pVec = this.v(1:end);
            else
                error('Empty queue.')
            end
        end
        
        function b = isEmpty(this)
            b = (this.n == 0);
        end
    end
end