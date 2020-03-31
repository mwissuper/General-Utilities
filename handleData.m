function handlePose = handleData(Markers,idx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%   Luke Drnach
%   Georgia Institue of Technology
%   March 26, 2018
    
%a row of NaNs causes the line between two points to NOT be drawn
skip = nan(1,3);

if isfield(Markers,{'frontleft','frontmiddle','frontright','backmiddle'})
   FL = Markers.frontleft(idx,:);
   if norm(FL) < 1
       FL = skip;
   end
   FM = Markers.frontmiddle(idx,:);
   if norm(FM) < 1
       FM = skip;
   end
   FR = Markers.frontright(idx,:);
   if norm(FR) < 1
       FR = skip;
   end
   BM = Markers.backmiddle(idx,:);
   if norm(BM) < 1
       BM = skip;
   end
   handlePose = [FL; FM; FR; BM; FL; FR; skip; FM; BM];
    
else
    handlePose = skip;
end

end

