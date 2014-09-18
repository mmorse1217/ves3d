precision = {'single','double'}
for p = [6,8,12,16] %[8 12 16 24 32 48 64 96 128]
  for ii=1:length(precision)
    disp(p);
    Generate_Cpp_data(p, precision{ii}, 'LegendreTransform', 'IntegrationWeights');
  
    if ( p <= 64 )
      Generate_Cpp_data(p, precision{ii}, 'SpHarmRot', 'DumbbellShape', ...
                        'SingularIntegralWeights', 'WSpherical','DirectRotation');
    end
  end
end