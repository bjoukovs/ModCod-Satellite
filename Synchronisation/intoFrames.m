function output = intoFrames(symb,frame_length, pilot)

    %Inserts a given pilot in symb each frame_length, at the beginning of
    %each frame

    pilot_length = length(pilot);
    symb_length = length(symb);

    nframes = symb_length/frame_length; %attention has to be an INTEGER, else the reshape will not work
                                        % because it will need a different
                                        % nbre of symbols than in symb
    frames = reshape(symb, [frame_length,nframes]);
    
    pilot = repmat(pilot, [1,nframes]);
    frames = [pilot; frames];
    
    output = reshape(frames,[symb_length + nframes*pilot_length, 1]);

end