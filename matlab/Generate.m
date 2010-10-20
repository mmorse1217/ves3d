precision = 'single';
for p=16
  disp(p);
  Generate_Cpp_data(p, precision, 'LegendreTransform', 'IntegrationWeights');
  
  if ( p <= 64 )
    Generate_Cpp_data(p, precision, 'SpHarmRot', 'DumbbellShape', ...
                      'SingularIntegralWeights', 'WSpherical','DirectRotation');
  end
end