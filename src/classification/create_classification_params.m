function params = create_classification_params( varargin )
    % params = create_classification_params( varargin )
    %
    % creates parameters for data classification. All inputs are optional,
    % and will default to the indicated value if not provided. To change,
    % please supply as name-value pairs
    %
    % inputs:
    %   fs - the sampling rate (default = 512 Hz)
    %   window - the amount (in seconds) of data to extract for each small
    %            chunk during classification (default = 0.25)
    %   wtype - the wavelet type (default = 'sym4')
    %   wlev - the wavelet level (default = 8)
    %   useSQRT - flag to sqrt transform the wavelet variances (default =
    %             false)
    %   useInteractions - flag to expand the wavelet features via
    %                     interactions (default = false)
    %
    % returns:
    %   params - a structure with the desired parameters
    
    
    params = input_parser;

    names = {'fs','window','wtype','wlev','useSQRT','useInteractions'};
    defaults = {512,0.25,'sym4',8,false,false};

    for j = 1:numel( names )
        params.addParameter(names{j},defaults{j});
    end

    parse( params,varargin{:} );
    params = params.Results;
end