classdef ContactStateMachine  < handle   
    properties (SetAccess = private, GetAccess = private)
        previousState (1, 1) logical;
        currentState (1, 1) logical;       
    end
    
    properties
        contactSchmitt Estimation.ContactHandler.SchmittTrigger;
    end
    
    methods
        function obj = ContactStateMachine(s)
            arguments                
                s Estimation.ContactHandler.SchmittParams;
            end
            
            obj = obj.setup(s);
        end
        
        function obj = setup(obj, s)
            obj.previousState = false;
            obj.currentState = false;
            
            obj.contactSchmitt = Estimation.ContactHandler.SchmittTrigger(s.stableTimeContactBreak, ...
                                                                          s.stableTimeContactMake, ...
                                                                          s.contactBreakForceThreshold, ...
                                                                          s.contactMakeForceThreshold);
            obj.contactSchmitt.setInitialState(obj.previousState);
        end
        
        function obj = contactMeasurementUpdate(obj, contactNormalForce, currentTime)
            arguments
                obj Estimation.ContactHandler.ContactStateMachine;
                contactNormalForce (1, 1) double;
                currentTime (1, 1) double;                
            end
            obj.contactSchmitt = obj.contactSchmitt.updateDevice(contactNormalForce, currentTime);
            obj.previousState = obj.currentState;
            obj.currentState = obj.contactSchmitt.getState();
        end
        
        function obj = resetDevice(obj)
            obj.contactSchmitt.resetDevice();
        end
        
        function contact_state = contactState(obj)
            contact_state = obj.currentState;
        end
        
        
        function contact_transition = contactTransitionMode(obj)
            if (obj.previousState == 0 && obj.currentState == 0)
                contact_transition = Estimation.ContactHandler.ContactTransition.STABLE_OFFCONTACT;                
                return
            end            
            if (obj.previousState == 0 && obj.currentState == 1)
                contact_transition = Estimation.ContactHandler.ContactTransition.CONTACT_MAKE;                
                return
            end            
            if (obj.previousState == 1 && obj.currentState == 0)
                contact_transition = Estimation.ContactHandler.ContactTransition.CONTACT_BREAK;                
                return
            end            
            if (obj.previousState == 1 && obj.currentState == 1)
                contact_transition = Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT;                
                return
            end                
            contact_transition = Estimation.ContactHandler.ContactTransition.UNKNOWN_TRANSITION;
        end
        
        function last_update_time = lastUpdateTime(obj)
            last_update_time = obj.contactSchmitt.getElapsedTime();
        end
        
    end
end
