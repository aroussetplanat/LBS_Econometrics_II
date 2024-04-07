function y=inv2(x)
row = size(x,1);
col = size(x,2);
    if ~isequal(row,col)
        errordlg('the matrix must be sequared')
    else
        newmat = [x,eye(row,col)];
%             i = 1;
        for i = 1:col %pivots
            pivot = newmat(i,i);
            yidash = i+1;
            if yidash > col
                yidash= yidash - col;
            end
            for j = i+1:i+col % col of pivot matrix
                %row of the pivot matrix
                rowyidash = yidash+ (row-2);
                if rowyidash > row
                    rowyidash = rowyidash - row;
                    yidashmat = [yidash:row 1:rowyidash];
                else
                    yidashmat = yidash:rowyidash;
                end
                
                for k = yidashmat % indices of pivot rows
                    newmat(k,j) = newmat(k,j)-((newmat(k,i)*newmat(i,j))/pivot);
                end
            end
            newmat(i,:)=newmat(i,:)/pivot;  
            newmat(:,i)=0;%pivot column
            newmat(i,i)=1;   
        end
    end
    y = newmat(:,col+1:end);
end