function R = MarginalProbability(bins_event_tble,events,numbins)

disp ('calculating marginal probabilities...');

%R=zeros(numbins,num_annot);
R=zeros(numbins,1);
for c1=1:length(events),
%     R(bins_event_tble(c1,1),(events(c1,3)-1)*2+events(c1,6))=R(bins_event_tble(c1,1),(events(c1,3)-1)*2+events(c1,6))+1;
%     R(bins_event_tble(c1,2),(events(c1,3)-1)*2+events(c1,6))=R(bins_event_tble(c1,2),(events(c1,3)-1)*2+events(c1,6))+1;
    R(bins_event_tble(c1,1))=R(bins_event_tble(c1,1))+1;
    R(bins_event_tble(c1,2))=R(bins_event_tble(c1,2))+1;
end

R = bsxfun(@rdivide,2*R,sum(R,1));


