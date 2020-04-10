% Given Array(1:n) returns the array Indx(1:n) such that Array(Indx(j)), j=1:n is in ascending order.
function Indx = indexArray(n, Array)

    Indx = zeros(1, n);
    
    sorted_Array = sort(Array);
    
    for i = 1 : length(Array)
        Indx(i) = find(Array == sorted_Array(i));
    end
    
end