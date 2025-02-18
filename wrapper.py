# Imports
from abc import ABC, abstractmethod
import subprocess
import threading
import os
import time

# Common
class IStartable():
    @abstractmethod
    def start(self) -> None:
        pass

class IStoppable():
    @abstractmethod
    def stop(self) -> None:
        pass

class IPipelineModelListener():
    @abstractmethod
    def onStateUpdate(self, state_id: int) -> None:
        pass

class IPipelineModelSource():
    @abstractmethod
    def setModelListener(listener: IPipelineModelListener) -> None:
        pass

class IPipelineSubprocessListener():
    @abstractmethod
    def on_start() -> None:
        pass
    def on_complete() -> None:
        pass

class ConcreteTimerTask():
    @abstractmethod
    def run() -> None:
        pass


class ConcreteTimer():
    def __init__(self):
        self.tasks = []
    
    def schedule(self, delay: int, interval: int, task: ConcreteTimerTask):
        def wrapper():
            task.run()
            timer = threading.Timer(interval, wrapper)
            self.tasks.append(timer)
            timer.start()
    
        first_timer = threading.Timer(delay, wrapper)
        self.tasks.append(first_timer)
        first_timer.start()
    
    def cancel(self):
        for timer in self.tasks:
            timer.cancel()
        self.tasks.clear()

# Model
## Subprocess
class IPollerListener:
    @abstractmethod
    def on_poll() -> None:
        pass

class IPollerSource:
    @abstractmethod
    def set_poller_listener(listener: IPollerListener) -> None:
        pass
        
class ISubprocessModel(IStartable, IStoppable, IPollerSource):
    @abstractmethod
    def _poll(self):
        pass

    @abstractmethod
    def _get_elapsed_time(self):
        pass

class DefaultSubprocessModel(ISubprocessModel):
    def __init__(self) -> None:
        self._listener: IPollerListener = IPollerListener()
        self._subprocess: subprocess.Popen | None = None
        self._running = False
        self._timer: threading.Timer | None = None
        self._start_time = 0
        self._current_time = 0

    def set_poll_listener(self, listener: IPollerListener) -> None:
        self._listener = listener    
    
    def _poll(self) -> None:
        if self._running and self._subprocess:
            retcode = self._subprocess.poll()
            if retcode is None:
                self._current_time = time.time()
                self._listener.on_poll()
            else:
                self._current_time = time.time()
                self._listener.on_complete()
                return   
        self._timer.start()
        
    def start(self, command: str) -> None:
        self.running = True
        self._subprocess = subprocess.Popen(os.sytem(command))
        timer = ConcreteTimer()
        task = ConcreteTimerTask()
        self._start_time = timer.time()
        task.run = lambda : self._listener.on_poll()
        timer.schedule(0, 5, task)
    
    def stop(self):
        self._running = False
        if self._timer:
            self._timer.cancel()
        if self._subprocess:
            self._subprocess.terminate()
            try:
                self._subprocess.wait(timeout = 2)
            except subprocess.TimeoutExpired:
                self._subprocess.kill()
            self.process = None
        self._current_time = time.time()

    
    def _get_elapsed_time(self) -> int:
        return int(self._current_time - self._start_time)
        

## State
class IPipelineState(IPipelineSubprocessListener, IPollerListener):
    @abstractmethod
    def _update_view() -> None:
        pass
    
    @abstractmethod
    def _get_id() -> int:
        pass

class IPipelineSMStateView:
    # Interactions with the subprocess model. 
    @abstractmethod
    def get_elapsed_time() -> int: pass

    # Transitions to other states. 
    @abstractmethod
    def to_mkdir_state() -> None: pass
    @abstractmethod
    def to_obtain_rnaseq_state() -> None: pass
    @abstractmethod 
    def to_trx_index_state() -> None: pass
    @abstractmethod
    def to_kallisto_state() -> None: pass
    @abstractmethod
    def to_sleuth_state() -> None: pass
    @abstractmethod
    def to_bowtie2_state() -> None: pass
    @abstractmethod
    def to_blast_state() -> None: pass

    # Actions.
    @abstractmethod
    def action_init() -> None: pass
    @abstractmethod
    def action_start() -> None: pass
    @abstractmethod
    def action_stop() -> None: pass
    @abstractmethod
    def write_to_log() -> None: pass
    @abstractmethod
    def print_progress() -> None: pass

class IPipelineStateMachine(IPipelineSubprocessListener, IPollerListener, IPipelineModelSource, IPipelineSMStateView):
    pass

class DefaultPipelineStateMachine(IPipelineStateMachine):
    def __init__(self, subprocess_model: ISubprocessModel):
        self._subprocess_model = subprocess_model
        self._state: IPipelineState | None = None
        self._listener: IPipelineModelListener | None = None
    
    def set_state(self, state: IPipelineState):
        self._state = state
        self._listener.onStateUpdate(state.get_id())
    
    def set_model_listener(self, listener: IPipelineModelListener):
        self._listener = listener

    def on_start(self):
        self._state.on_start()
    
    def on_complete(self):
        self._state.on_complete()


    
