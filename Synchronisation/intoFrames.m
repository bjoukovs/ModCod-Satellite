function output = intoFrames(symb,frame_length, pilot)

    %Inserts a given pilot in symb each frame_length, at the beginning of
    %each frame

    pilot_length = length(pilot);
    symb_length = length(symb);

    nframes = ceil(symb_length/frame_length);
    frames = reshape(symb, [frame_length,nframes]);
    
    pilot = repmat(pilot, [1,nframes]);
    frames = [pilot; frames];
    
    output = reshape(frames,[symb_length + nframes*pilot_length, 1]);

end