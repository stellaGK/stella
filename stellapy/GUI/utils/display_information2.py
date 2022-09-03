def display_information(root,title=None,message=None):
    ''' Display a window with information.

    Parameters
    ----------
    self : class object
        self.root should be linked to the root window in order to center the window on the screen.
    message : str
        Message to be displayed.
    '''

    if not title:   title="Please insert a title"
    if not message: message="Please insert a message."
    information_window = tk.Toplevel(bg=root.color['bg'])
    information_window.title(title)
    root.eval(f'tk::PlaceWindow {str(information_window)} center')
    return
