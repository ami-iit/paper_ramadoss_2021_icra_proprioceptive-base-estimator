classdef ContactTransition
    enumeration
        STABLE_OFFCONTACT, ...
        STABLE_ONCONTACT, ...
        CONTACT_BREAK, ...
        CONTACT_MAKE, ...
        UNKNOWN_TRANSITION
    end
    
    methods (Static)
        function contact_transition = contactTransitionMode(previousState, currentState)
            if (previousState == 0 && currentState == 0)
                contact_transition = Estimation.ContactHandler.ContactTransition.STABLE_OFFCONTACT;
                return
            end
            if (previousState == 0 && currentState == 1)
                contact_transition = Estimation.ContactHandler.ContactTransition.CONTACT_MAKE;
                return
            end
            if (previousState == 1 && currentState == 0)
                contact_transition = Estimation.ContactHandler.ContactTransition.CONTACT_BREAK;
                return
            end
            if (previousState == 1 && currentState == 1)
                contact_transition = Estimation.ContactHandler.ContactTransition.STABLE_ONCONTACT;
                return
            end
            contact_transition = Estimation.ContactHandler.ContactTransition.UNKNOWN_TRANSITION;
        end
    end
end

