classdef ActiveSet < handle
    properties
        tol = 1e-8;
        maxiter = 10000;
        neq = 0;
        nbds = 0;
        A 
        b 
        C
        d 
        AC_KKT
        acidx 
        H 
        h 
        Ce
        de 
        ce
        cl
        cu
    end

    methods
        function obj = ActiveSet(atol, maxiter)
            if nargin >= 1
                obj.tol = atol;
            end
            if nargin >= 2
                obj.maxiter = maxiter;
            end
        end

        function init_vars(obj)
            obj.neq = 0;
            obj.nbds = 0;
            obj.A = [];
            obj.b = [];
            obj.C = [];
            obj.d = [];
            obj.AC_KKT = [];
            obj.acidx = [];
            obj.H = [];
            obj.h = [];
            obj.Ce = [];
            obj.ce = [];
            obj.cl = [];
            obj.cu = [];
        end

        function calc_Hessians(obj)
            obj.H = obj.A;
            obj.h = obj.b;
        end

        function hvec = calc_as_Hessian_vector(obj, x)
            hvec = obj.h - obj.H*x;
        end

        function val = calc_objective(obj, x)
            val = 0.5*x'*obj.A*x - obj.b'*x;
        end

        function obj = add_active_constraint(obj, cidx)
            % Add the constraint index to acidx
            obj.acidx = [obj.acidx; cidx];

            % Extract the constraint row vector
            c_row = obj.C(cidx, :);

            % Initialize new row to append to AC_KKT
            new_row = zeros(1, size(obj.AC_KKT, 2));
            new_row(1:length(c_row)) = c_row;

            % Append the new row to AC_KKT
            obj.AC_KKT = [obj.AC_KKT; new_row];

            % Initialize new column to append to AC_KKT
            new_col = zeros(size(obj.AC_KKT, 1), 1);
            new_col(1:length(c_row)) = c_row';

            % Append the new column to AC_KKT
            obj.AC_KKT = [obj.AC_KKT, new_col];
        end


        function obj = remove_active_constraint(obj, wcidx)
            % Remove the constraint index from acidx
            obj.acidx(wcidx) = [];

            % Calculate the KKT index to remove (offset by size of H)
            kktidx = wcidx + size(obj.H, 1);

            % Remove the corresponding row and column from AC_KKT
            obj.AC_KKT(kktidx, :) = [];
            obj.AC_KKT(:, kktidx) = [];
        end


        function [alpha, cidx] = compute_step_length(obj, x, p)
            % constraint indices
            cidx_all = (1:size(obj.C, 1))';
            
            % get inactive constraints
            mask = ~ismember(cidx_all, obj.acidx);
            IC = obj.C(mask, :);
            id = obj.d(mask);
            
            % keep track of constraint matrix indices
            icidx = cidx_all(mask);
            
            % denominator
            den = IC * p;
            
            % avoid division by zero
            q = (id - IC * x) ./ den;
            
            % filter for positive denominators
            valid = den > obj.tol;
            q = q(valid);
            icidx = icidx(valid);
            
            % check if empty
            if isempty(q)
                alpha = 1;
                cidx = -1;
                return;
            end
            
            [alpha, minidx] = min(q);
            cidx = icidx(minidx);
        end




        function [p, l] = null_space_as_kkt(obj, x)
            % Number of active constraints
            m = size(obj.AC_KKT, 1) - size(obj.H, 1);

            % If there are no active constraints
            if m == 0
                p = obj.H \ obj.calc_as_Hessian_vector(x);
                l = [];
                return;
            end

            % Number of variables
            n = size(obj.C, 2);

            % Zero right-hand side for active constraints
            dim = zeros(m, 1);

            % Compute Hessian-vector product
            hessi = obj.calc_as_Hessian_vector(x);

            % QR decomposition of the active constraints (transpose)
            [Q, R_full] = qr(obj.AC_KKT(1:size(obj.H,1), size(obj.H,2)+1:end), 0);

            % Trim R to upper m x m
            R = R_full(1:m, 1:m);

            % Solve R * py = d
            py = R \ dim;

            % Compute range-space projection
            Qr = Q(:, 1:m);
            p = Qr * py;

            % If we are in the null space, compute reduced step
            if n > m
                Qn = Q(:, m+1:end);
                Qnt = Qn';

                % Reduced Hessian
                Hz = Qnt * obj.H * Qn;

                % Cholesky factorization
                L = chol(Hz, 'lower');

                % Compute null-space target vector
                hz = Qnt * (obj.H * Qr * py + hessi);

                % Solve for null-space step
                pz = L' \ (L \ hz);

                % Add null-space component to step
                p = p + Qn * pz;
            end

            % Solve for Lagrange multipliers
            l = R \ (Qr' * (obj.H * (-p) + hessi));
        end


        function [p, l] = solve_as_kkt(obj, x)
            % Check if active constraint set is empty
            if size(obj.AC_KKT, 1) - size(obj.H, 1) == 0
                p = obj.H \ obj.calc_as_Hessian_vector(x);
                l = [];
                return;
            end
        
            % Construct the KKT target vector using current solution
            kkt = zeros(size(obj.AC_KKT, 1), 1);
            kkt(1:size(obj.A, 2)) = obj.calc_as_Hessian_vector(x);
        
            % Solve the KKT system
            I = eye(size(obj.AC_KKT));
            wx_wl = (obj.AC_KKT+ eps*I) \ kkt;
        
            % Extract step direction and Lagrange multipliers
            p = wx_wl(1:size(obj.A, 2));
            l = wx_wl(size(obj.A, 2)+1:end);
        end


        function init_active_set(obj, x0)
            % Initialize the KKT matrix with the Hessian
            obj.AC_KKT = obj.H;
        
            % Compute current constraint values
            cx = obj.C * x0;
        
            % Initialize the active constraint index list
            obj.acidx = zeros(0, 1);
        
            % Loop through all constraints
            for j = 1:length(cx)
                if j <= obj.neq || abs(cx(j) - obj.d(j)) <= obj.tol
                    obj.add_active_constraint(j);
                end
            end
        end


        function check_1d_array(obj, m, mnm)
            if ndims(m) > 2 || (ismatrix(m) && size(m, 2) ~= 1)
                error('ActiveSet: %s must be a 1D array or 2D with a single column', mnm);
            end
        end


        function [M_out, m_out] = check_constraints(obj, M, m, mtxnm, mnm, mtxemptyok)
            if nargin < 6
                mtxemptyok = true;
            end
        
            M = double(M);
            m = double(m);
        
            % Objective matrix cannot be empty
            if ~mtxemptyok && isempty(M)
                error('ActiveSet: %s is empty', mtxnm);
            elseif ~isempty(M)
                % Matrix must be 2D
                if ndims(M) ~= 2
                    error('ActiveSet: %s must be a 2D array', mtxnm);
                end
                % Constraint matrices must match A column size
                if ~strcmp(mtxnm, 'A') && size(M, 2) ~= size(obj.A, 2)
                    error('ActiveSet: %s (axis=1) must match A (axis=1)', mtxnm);
                end
            end
        
            % Check m is 1D or a column vector
            obj.check_1d_array(m, mnm);
        
            % Ensure matrix-vector pairs have matching dimensions
            if strcmp(mnm, 'b') && ~isempty(m) && size(M, 1) ~= size(m, 1)
                error('ActiveSet: b must be empty or match A (axis=0)');
            elseif ~strcmp(mnm, 'b') && size(M, 1) ~= size(m, 1)
                error('ActiveSet: %s and %s mismatch (axis=0)', mtxnm, mnm);
            end
        
            % Reshape if empty
            if isempty(M)
                M = zeros(0, size(obj.A, 2));
            end
            m = reshape(m, [], 1);
        
            % Append to C and d if not for A
            if ~strcmp(mtxnm, 'A') && ~isempty(M)
                obj.C = [obj.C; M];
                obj.d = [obj.d; m];
            end
        
            M_out = M;
            m_out = m;
        end


        function bds_out = check_bounds(obj, b, bnm, upper)
            if nargin < 4
                upper = true;
            end
        
            bds = double(b);
        
            if ~isempty(bds)
                obj.check_1d_array(bds, bnm);
                if size(bds, 1) ~= size(obj.A, 2)
                    error('ActiveSet: bounds (axis=0) must match A (axis=1)');
                end
            end
        
            bds = reshape(bds, [], 1);
        
            if upper && ~isempty(bds)
                obj.C = [obj.C; eye(size(bds, 1))];
                obj.d = [obj.d; bds];
            elseif ~upper && ~isempty(bds)
                obj.C = [obj.C; -eye(size(bds, 1))];
                obj.d = [obj.d; -bds];
            end
        
            bds_out = bds;
        end

        function prep_inputs(obj, A, b, Ce, de, cu, cl)
            % initialize instance variables
            obj.init_vars();
        
            % set objective function instance variables
            % objective matrix cannot be empty. target vector can be empty
            [obj.A, obj.b] = obj.check_constraints(A, b, 'A', 'b', false);
        
            % initialize C, d
            obj.C = reshape(obj.C, [], size(obj.A, 2));
            obj.d = reshape(obj.d, [], 1);
        
            % set Hessian instance variables
            obj.calc_Hessians();
        
            % add equality constraints
            [obj.Ce, obj.de] = obj.check_constraints(Ce, de, 'Ce', 'de');
            obj.neq = size(obj.Ce, 1);
        
            % add upper bound constraints
            obj.cu = obj.check_bounds(cu, 'cu');
            obj.nbds = size(obj.cu, 1);
        
            % add lower bound constraints
            obj.cl = obj.check_bounds(cl, 'cl', false);
            obj.nbds = obj.nbds + size(obj.cl, 1);
        
            % flag empty input
            if isempty(obj.C)
                error(['ActiveSet: input must have at least one of ', ...
                       '[eq,ineq,bound] constraints']);
            end
        end


        function [cur_x, fval, iter] = run(obj, A, b, Ce, de, cu, cl, x0)
            % Construct input
            obj.prep_inputs(A, b, Ce, de, cu, cl);
        
            obj.init_active_set(x0);
        
            % Main loop
            cur_x = x0;
            for iter = 1:obj.maxiter
                % Solve for step direction and Lagrangian multipliers
                try
                    [p, l] = obj.solve_as_kkt(cur_x);
                catch
                    [p, l] = obj.null_space_as_kkt(cur_x);
                end
        
                % Check length of step
                len_p = norm(p);
                if abs(len_p - 0) < obj.tol
                    if size(l, 1) == obj.neq
                        fval = obj.calc_objective(cur_x);
                        return;
                    end
        
                    m = min(l(obj.neq+1:end));
                    if m >= -obj.tol
                        fval = obj.calc_objective(cur_x);
                        return;
                    end
        
                    nlidx = find(l == m, 1);
                    obj.remove_active_constraint(nlidx);
                else
                    % Calculate step length
                    [alpha, cidx] = obj.compute_step_length(cur_x, p);
        
                    if abs(alpha - 0) <= obj.tol
                        fval = obj.calc_objective(cur_x);
                        return;
                    elseif alpha < 1
                        obj.add_active_constraint(cidx);
                        cur_x = cur_x + alpha * p;
                    else
                        cur_x = cur_x + p;
                    end
                end
            end
        
            % Issue max iteration warning
            warning('ActiveSet: maximum number of iterations reached (%d)', obj.maxiter);
            fval = obj.calc_objective(cur_x);
        end

        function result = isclose(a, b, tol)
            if nargin < 3
                tol = 1e-8;
            end
            result = abs(a - b) <= tol;
        end

    end
end


