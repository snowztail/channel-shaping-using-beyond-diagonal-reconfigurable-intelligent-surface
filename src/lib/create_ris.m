function [ris] = create_ris(nSxs, nGroups)
	ris.scatter = eye(nSxs);
    ris.group = nGroups;
    ris.connect = nSxs / nGroups;
end
