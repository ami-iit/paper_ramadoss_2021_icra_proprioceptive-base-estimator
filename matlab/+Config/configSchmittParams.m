function [leftFootSchmittParams, rightFootSchmittParams] = configSchmittParams()
leftFootSchmittParams = Estimation.ContactHandler.SchmittParams();
rightFootSchmittParams = Estimation.ContactHandler.SchmittParams();

leftFootSchmittParams.contactBreakForceThreshold = 120.0;
rightFootSchmittParams.contactBreakForceThreshold = 120.0;

leftFootSchmittParams.contactMakeForceThreshold = 150.0;
rightFootSchmittParams.contactMakeForceThreshold = 150.0;

leftFootSchmittParams.stableTimeContactBreak = 0.01; 
rightFootSchmittParams.stableTimeContactBreak = 0.01; 

leftFootSchmittParams.stableTimeContactMake = 0.01; 
rightFootSchmittParams.stableTimeContactMake = 0.01; 


end

