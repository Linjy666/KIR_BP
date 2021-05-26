classdef Player < handle
    properties
        state
    end
    properties(Constant)
        PLAYER_STOP = 0;
        PLAYER_PLAYING = 1;
        PLAYER_PAUSE = 2;        
    end
    methods
        function obj = Player(state)
            obj.state = state;
        end
        function request(obj, state)
            obj.state.handleRequest(obj);
        end
%         function stop(obj, state)
%             obj.state.stop(obj);
%         end
%         function pause(obj, state)
%             obj.state.pause(obj);
%         end
    end
end