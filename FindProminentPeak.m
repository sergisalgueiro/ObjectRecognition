function idxs = FindProminentPeak(histogram, width)

peaks = [];

    for k=width+1:length(histogram)-width

        pattern=histogram([k-width k k+width])';

        [~,ind]=max(histogram(k-width:k+width));
        ind=ind(1);

        %val = pattern(2) - 2*pattern(1) - 2*pattern(end);
        val = pattern(2) - pattern(1) - pattern(end);

        if((ind == width+1) && k < (length(histogram)))
            peaks = [peaks; [k val]];
        end
    end

idxs = peaks;

end