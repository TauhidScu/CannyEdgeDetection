% Find weak edges that are connected to strong edges and set them to 1
function[Gmag] = FindConnectedWeakEdges(Gmag, row, col)
    for i = -3:1:3
        for j = -3:1:3
            if (row+i > 0) && (col+j > 0) && (row+i < size(Gmag,1)) && ...
                    (col+j < size(Gmag,2)) % Make sure we are not out of bounds
                if (Gmag(row+i,col+j) > 0) && (Gmag(row+i,col+j) < 1)
                    Gmag(row+i,col+j) = 1;
                    Gmag = FindConnectedWeakEdges(Gmag, row+i, col+j);
                end
            end
        end
    end
end

