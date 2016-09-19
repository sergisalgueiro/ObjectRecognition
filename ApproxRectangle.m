function output = ApproxRectangle(input)
% It determines the 4 corners of the tree approximation.
% Sergi Salgueiro
output =[min(input(:,1)) min(input(:,2));
         min(input(:,1)) max(input(:,2));
         max(input(:,1)) max(input(:,2));
         max(input(:,1)) min(input(:,2))];
end