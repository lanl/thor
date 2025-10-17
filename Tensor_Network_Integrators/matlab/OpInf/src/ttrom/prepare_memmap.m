function mmap = prepare_memmap(cfg)
  % total dimensions
  mmap.Nt    = cfg.Nt;
  mmap.rows  = cfg.N * cfg.Neq;

  % how many bytes per full column?
  bytesPerCol = mmap.rows * 8;         

  % maximum number of columns that fit in targetGB
  mmap.chunkCols = max(1, floor(cfg.targetGB * 2^30 / bytesPerCol));
  mmap.nChunks   = ceil(mmap.Nt / mmap.chunkCols);

  % preallocate metadata
  mmap.startCol = zeros(mmap.nChunks,1);
  mmap.colCount = zeros(mmap.nChunks,1);
  mmap.maps     = cell(mmap.nChunks,1);

  % build each time-chunk
  for cc = 1:mmap.nChunks
    % which columns this chunk covers
    mmap.startCol(cc) = (cc-1)*mmap.chunkCols + 1;
    thisCols          = min(mmap.chunkCols, mmap.Nt - (cc-1)*mmap.chunkCols);
    mmap.colCount(cc) = thisCols;

    % byte‐offset = header + (# full columns before) * bytesPerCol
    offset = cfg.offset + (cc-1)*bytesPerCol*mmap.chunkCols;

    % map [rows × thisCols] doubles in one contiguous block
    mmap.maps{cc} = memmapfile(cfg.filename, ...
      'Offset',   offset, ...
      'Format',   {'double',[mmap.rows, thisCols],'Q'}, ...
      'Repeat', 1, ...
      'Writable', false);
  end
end
